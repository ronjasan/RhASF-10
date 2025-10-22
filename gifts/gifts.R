#### Script to make plots plot Figure 2

# Author: Ronja Marlonsdotter Sandholm

library(tidyverse)
library(readxl)
library(distillR)
library(Rtsne)
library(Polychrome)
library(broom)
library(ggstar)
library(ggrepel)
library(cowplot)
library(pheatmap)

## GIFT calculation
annotations_control <- read_tsv("data/cantarcticum_annotations.tsv") %>%
    mutate(Accession = ifelse(grepl("GCA", fasta), str_extract(fasta, "GCA_[0-9]+\\.[0-9]+"), "ASF-10")) %>%
    mutate(Gene_id = paste0(scaffold, "_", gene_position)) %>%
    select(Accession, Gene_id, ko_id)

annotations <- read_tsv("data/Rhodococcus_annotations.tsv") %>%
    mutate(Accession = ifelse(grepl("GCA", fasta), str_extract(fasta, "GCA_[0-9]+\\.[0-9]+"), "ASF-10")) %>%
    mutate(Gene_id = paste0(scaffold, "_", gene_position)) %>%
    select(Accession, Gene_id, ko_id) %>%
    bind_rows(annotations_control)

GIFTs <- distill(annotations, GIFT_db, genomecol = 1, annotcol = 3)

GIFTs_elements <- to.elements(GIFTs, GIFT_db)
GIFTs_functions <- to.functions(GIFTs_elements, GIFT_db)
GIFTs_domains <- to.domains(GIFTs_functions, GIFT_db)

avg_mci <- rowMeans(GIFTs_functions) %>%
    enframe(name = "Accession", value = "MCI")

metadata <- read_xlsx("../tree/data/strain_info.xlsx", sheet = 1) %>%
    mutate(Accession = ifelse(is.na(Accession), "ASF-10", Accession)) %>%
    mutate(Clade = as.factor(Clade), Plastic = fct_relevel(Plastic, "Polyethylene", "Crude oil", "Alkane", "No")) %>%
    mutate(
        Source = fct_relevel(Source, "Host: animal", "Host: plant", "Soil", "Water", "Other"),
        Species_short = str_replace(Species, "^Rhodococcus(?!\\s+sp\\.)", "R.")
    )

avg_mci_meta <- avg_mci %>%
    left_join(metadata, by = "Accession")

## t-SNE calculation
set.seed(123)

tsne_func <- Rtsne(
    X = GIFTs_elements,
    dims = 2,
    perplexity = min(30, floor((nrow(GIFTs_elements) - 1) / 3)),
    check_duplicates = FALSE
)

tsne_df <- as_tibble(tsne_func$Y) %>%
    rename(tSNE1 = V1, tSNE2 = V2) %>%
    bind_cols(avg_mci_meta) %>%
    select(Accession, everything())

TSNEplotOrder <- tsne_df %>%
    ggplot(aes(
        x = tSNE1,
        y = tSNE2,
        fill = Clade,
        shape = Plastic
    )) +
    # geom_point(size = 2, alpha = 0.8) +
    geom_star(aes(starshape = Plastic, fill = Clade), size = 3, color = "black") +
    scale_fill_manual(values = c("#ff7f00", "#fdae6b", "#7D3200", "#B78560", "#969696", "#d9d9d9", "#1c8457")) +
    geom_text_repel(aes(label = Species_short),
        size = 3, color = "black",
        show.legend = F,
        hjust = 1,
        vjust = 2
    ) +
    theme_classic() +
    theme(
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.75),
        axis.text = element_text(color = "black", size = 14),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(margin = margin(5, 0, 0, 0)),
        axis.text.y = element_text(margin = margin(0, 5, 0, 0)),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )


mci_color <- colorRampPalette(c("white", "#fdeda0", "#8452a1"))(10)
TSNEplotMCI <- tsne_df %>%
    ggplot(aes(
        x = tSNE1,
        y = tSNE2,
        color = MCI
    )) +
    geom_star(aes(starshape = Plastic, fill = MCI), size = 3, color = "black") +
    scale_fill_gradientn(
        colors = mci_color
    ) +
    geom_text_repel(aes(label = Species_short),
        size = 3, color = "black",
        show.legend = F,
        hjust = 1,
        vjust = 2
    ) +
    theme_classic() +
    theme_classic() +
    theme(
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.75),
        axis.text = element_text(color = "black", size = 14),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(margin = margin(5, 0, 0, 0)),
        axis.text.y = element_text(margin = margin(0, 5, 0, 0)),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )

plot <- plot_grid(
    TSNEplotOrder,
    TSNEplotMCI,
    ncol = 1,
    labels = "AUTO",
    label_size = 20
)

# Heatmap of GIFTs
Elements <- GIFT_db %>%
    select(Code_element, Domain, Function)

GIFTs_elements %>%
    pheatmap::pheatmap(cluster_cols = F)

GIFTs_W_Tax <- as_tibble(GIFTs_elements, rownames = NA) %>%
    rownames_to_column("Accession") %>%
    left_join(Metadata, by = "Accession")

GIFTSannotCol <- Elements %>%
    select(Code_element, Function) %>%
    distinct() %>%
    column_to_rownames("Code_element")


GIFTSTaxAnnot <- GIFTs_W_Tax %>%
    select(Accession, Clade) %>%
    column_to_rownames("Accession")

GIFTSPlasticAnnot <- GIFTs_W_Tax %>%
    select(Accession, Plastic) %>%
    column_to_rownames("Accession")


ColorAnnot <- list(
    Clade = c(
        "1" = "#ff7f00",
        "2" = "#fdae6b",
        "3" = "#7D3200",
        "4" = "#B78560",
        "5" = "#969696",
        "6" = "#d9d9d9",
        "Control" = "#1c8457"
    ),
    Plastic = c(
        "Polyethylene" = "#ff7f00",
        "Crude oil" = "#7D3200",
        "Alkane" = "#969696",
        "No" = "#000000"
    ),
    Function = (Elements %>%
        select(Function) %>%
        distinct() %>%
        mutate(Color = palette36.colors(n = 21)) %>%
        deframe())
)

color_heat <- colorRampPalette(c("#feeda0", "#f57c00", "#7D3200"))(10)

GIFTsElementsPH <- GIFTs_elements %>%
    pheatmap::pheatmap(
        cluster_cols = F,
        color = color_heat,
        annotation_row = GIFTSTaxAnnot,
        annotation_col = GIFTSannotCol,
        annotation_colors = ColorAnnot,
        show_colnames = F,
        show_rownames = F,
        annotation_legend = FALSE
    )


# ANOVA
aov <- tsne_df %>%
    mutate(size = case_when(
        Genome_size < 6000000 ~ "<6 Mb",
        Genome_size >= 6000000 & Genome_size < 8000000 ~ "6-8 Mb",
        Genome_size > 8000000 ~ ">8 Mb"
    )) %>%
    select(Accession, MCI, size) %>%
    column_to_rownames("Accession")

mci_avg <- aov %>%
    group_by(size) %>%
    summarise(mean_MCI = mean(MCI), sd_MCI = sd(MCI), n = n())

## One-way ANOVA
anova_result <- aov(MCI ~ size, data = aov)
summary(anova_result)

# Post-hoc tests if ANOVA is significant
if (summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05) {
    # Tukey's HSD for pairwise comparisons
    tukey_result <- TukeyHSD(anova_result)
    print(tukey_result)
}
