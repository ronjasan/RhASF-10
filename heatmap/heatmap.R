#### Script to make heatmap Figure 2

# Author: Ronja Marlonsdotter Sandholm
library(tidyverse)
library(microshades)
library(readxl)
library(DT)

blast <- read_tsv("heatmap/data/blast_hits.txt", col_names = c("qseqid", "sseqid", "qlen", "length", "pident", "qcovhsp", "bitscore", "evalue")) %>%
    filter(pident >= 50, qcovhsp >= 75, evalue <= 1e-5) %>%
    mutate(
        accession = str_extract(sseqid, "^[^|]*"),
        qseqid = str_replace(qseqid, "2_01_Rhodococcus_sp002259485", "RhASF-10"),
        rhodo = case_when(
            grepl("RhASF-10", qseqid) ~ "RhASF-10",
            grepl("GCA", qseqid) ~ str_remove(qseqid, "\\|.*"),
            TRUE ~ qseqid
        )
    ) %>%
    mutate(id = qseqid) %>%
    mutate(
        qseqid = case_when(
            grepl("RhASF-10", qseqid) ~ paste(str_extract(qseqid, "^[^_]+"), sseqid, sep = "_"),
            grepl("GCA", qseqid) ~ paste(str_extract(qseqid, "^[^|]*"), sseqid, sep = "_"),
            TRUE ~ qseqid
        )
    ) %>%
    arrange(desc(pident)) %>%
    group_by(qseqid) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(desc(pident)) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()

info <- read_xlsx("tree/data/strain_info.xlsx") %>%
    rename(rhodo = Accession) %>%
    mutate(rhodo = ifelse(is.na(rhodo), "RhASF-10", rhodo))

meta_alkanes_plastics <- read_tsv("heatmap/data/250915_plastic_db_meta.tsv") %>%
    right_join(blast, by = "accession") %>%
    mutate(
        rhodo = ifelse(grepl("GCA_", qseqid), str_extract(qseqid, "^[^_]+_[^_]+"), "RhASF-10"),
        gene = case_when(
            gene == "lac" ~ "lacc",
            TRUE ~ gene
        )
    ) %>%
    mutate(key = paste(gene, accession, sep = "_")) %>%
    group_by(gene) %>%
    mutate(key_2 = paste0(gene, "_", dense_rank(key))) %>%
    ungroup() %>%
    mutate(pathway = case_when(
        compound == "Plastics" & grepl("/", pathway) ~ "Multiple",
        pathway == "LDPE" ~ "PE",
        pathway == "Nylon" ~ "PA",
        TRUE ~ pathway
    )) %>%
    mutate(pathway = fct_relevel(pathway, c("Auxiliar_alkane_gene", "Terminal/biterminal_oxidation", "Subterminal_oxidation", "Finnerty_pathway", "PE", "PET", "PUR", "PA", "PEG", "PHA", "NR", "Multiple"))) %>%
    full_join(info, by = "rhodo") %>%
    mutate(Species = fct_relevel(Species, c(
        "Corynebacterium antarcticum CCM 8835",
        "Rhodococcus ruber C1",
        "Rhodococcus xishaensis LHW51113",
        "Rhodococcus maanshanensis DSM 44675",
        "Rhodococcus spelaei C9-5",
        "Rhodococcus trifolii CCM 7905",
        "Rhodococcus sp. ASF-10",
        "Rhodococcus jostii DSM 44719",
        "Rhodococcus koreensis DSM 44498",
        "Rhodococcus opacus DSM 43205",
        "Rhodococcus opacus R7",
        "Rhodococcus wratislaviensis NBRC100605",
        "Rhodococcus globerulus NBRC14531",
        "Rhodococcus qingshengii RL1",
        "Rhodococcus qingshengii TG-1",
        "Rhodococcus qingshengii CL-05",
        "Rhodococcus qingshengii JCM15477",
        "Rhodococcus qingshengii BF1",
        "Rhodococcus qingshengii A34",
        "Rhodococcus sp. C-2",
        "Rhodococcus erythropolis R138",
        "Rhodococcus qingshengii 7B"
    )))

datatable(meta_alkanes_plastics)

plot <- ggplot(meta_alkanes_plastics, aes(x = key_2, y = Species)) +
    geom_tile(aes(fill = pident), color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
        high = microshades_palettes$micro_brown[5],
        mid = microshades_palettes$micro_orange[5],
        low = microshades_palettes$micro_orange[1],
        midpoint = 75,
        limits = c(50, 100)
    ) +
    facet_grid(cols = vars(compound, pathway), scales = "free_x", space = "free") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold"),
        # legend.position = "bottom",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.75)
    )
