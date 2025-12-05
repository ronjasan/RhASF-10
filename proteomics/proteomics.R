#### Script to make volcano plot Figure 5

# Author: Ronja Marlonsdotter Sandholm

library(tidyverse)
library(rstatix)
library(DT)
library(rstatix)
library(ggrepel)
library(cowplot)

## Plastic genes
blast <- read_tsv("data/plastic_db_blast.txt", col_names = c("qseqid", "sseqid", "qlen", "length", "pident", "qcovhsp", "bitscore", "evalue")) %>%
    filter(pident >= 50, qcovhsp >= 75, evalue <= 1e-5) %>%
    mutate(accession = str_remove(sseqid, "\\|.*")) %>%
    mutate(id = qseqid) %>%
    mutate(
        qseqid = case_when(
            grepl("RhRMS", qseqid) ~ paste(str_extract(qseqid, "^[^_]+"), sseqid, sep = "_"),
            grepl("GCA", qseqid) ~ paste(str_extract(qseqid, "^([^_]+_[^_]+)"), sseqid, sep = "_"),
            TRUE ~ qseqid
        )
    ) %>%
    arrange(desc(pident)) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(desc(pident)) %>%
    group_by(qseqid) %>%
    slice(1) %>%
    ungroup() %>%
    filter(grepl("2_01_Rhodococcus", qseqid))


meta <- read_tsv("data/plastic_db_meta.tsv") %>%
    distinct() %>%
    right_join(blast, by = "accession") %>%
    select(compound, pathway, gene, enzyme, pident, id)


# Annotations

DRAM <- read_tsv("data/RhASF-10_annotations.tsv") %>%
    rename(id = `...1`) %>%
    left_join(meta, by = "id") %>%
    mutate(id = str_replace(id, "2_01_Rhodococcus_sp002259485", "RhRMS")) %>%
    mutate(dram_annotation = case_when(rank == "E" ~ "unknown", rank == "D" ~ pfam_hits, rank == "C" ~ kegg_hit)) %>%
    mutate(dram_annotation = str_replace(dram_annotation, "(\\]; ).*", "]")) %>%
    mutate(pathway = case_when(
        compound == "Alkanes" ~ "alkane",
        compound == "Plastics" ~ "plastic",
        TRUE ~ NA_character_
    )) %>%
    select(id, KO = ko_id, gene1 = gene, pathway1 = pathway, enzyme, dram_annotation)


KEGG <- read_tsv("data/KEGG_pathways.tsv") %>%
    right_join(DRAM, by = "KO") %>%
    mutate(gene = case_when(
        is.na(gene) & !is.na(gene1) ~ gene1,
        gene == "fmo" & gene1 == "almA" ~ "almA",
        gene == "cyp" & gene1 == "CYP153" ~ "cyp153",
        is.na(gene) & is.na(KO) & grepl("Cytochrome P450", dram_annotation) ~ "cyp",
        is.na(gene) & is.na(KO) & grepl("Rubredoxin", dram_annotation) ~ "rub",
        is.na(gene) & is.na(KO) & grepl("Flavin-binding monooxygenase-like", dram_annotation) ~ "bvmo",
        is.na(gene) & is.na(KO) & grepl("Alcohol dehydrogenase GroES-like domain", dram_annotation) ~ "adh",
        is.na(gene) & is.na(KO) & grepl("Carboxylesterase family", dram_annotation) ~ "ces",
        TRUE ~ gene
    )) %>%
    mutate(pathway = case_when(
        is.na(pathway) & !is.na(pathway1) ~ pathway1,
        is.na(pathway) & grepl("cyp", gene) ~ "alkane",
        is.na(pathway) & grepl("rub", gene) ~ "alkane",
        is.na(pathway) & grepl("bvmo", gene) ~ "ketone",
        is.na(pathway) & grepl("ces", gene) ~ "ketone",
        is.na(pathway) & grepl("adh", gene) ~ "alkane",
        TRUE ~ pathway
    )) %>%
    select(id, gene, pathway, dram_annotation) %>%
    mutate(
        protein = case_when(
            str_detect(gene, "[A-Z0-9]") ~ paste0(toupper(substr(gene, 1, 1)), substr(gene, 2, nchar(gene))),
            TRUE ~ toupper(gene)
        )
    ) %>%
    mutate(protein = if_else(!is.na(protein), paste0(protein, "_", str_extract(id, "(?<=_).*")), protein))


# Proteomics
source("impute_normal.R")

intensity <- read_tsv("data/RhASF-10_combined_protein.tsv") %>%
    select(Protein, matches("^(Succinate|LMWPE)_[1-3] Intensity$")) %>%
    pivot_longer(cols = -c(Protein), names_to = "Sample", values_to = "Intensity") %>%
    mutate(Sample = str_remove(Sample, " Intensity")) %>%
    mutate(Intensity = log2(Intensity)) %>%
    mutate(Intensity = if_else(is.infinite(Intensity), NA_real_, Intensity)) %>%
    pivot_wider(names_from = Sample, values_from = Intensity) %>%
    filter(
        rowSums(is.na(select(., matches("^Succinate_[1-3]$")))) <= 1 &
            rowSums(is.na(select(., matches("^LMWPE_[1-3]$")))) <= 1
    ) %>%
    mutate(across(where(is.numeric), ~ impute_normal(.))) %>%
    mutate(across(where(is.numeric), ~ as.numeric(.))) %>%
    rename(id = Protein) %>%
    filter(!grepl("sp", id))

diff <- intensity %>%
    mutate(
        succ_mean = rowMeans(select(., matches("^Succinate_[1-3]$")), na.rm = TRUE),
        LMWPE_mean = rowMeans(select(., matches("^LMWPE_[1-3]$")), na.rm = TRUE),
        diff = LMWPE_mean - succ_mean
    )

ttest <- intensity %>%
    pivot_longer(cols = -id, names_to = "Sample", values_to = "Intensity") %>%
    mutate(
        condition = case_when(
            str_detect(Sample, "Succinate") ~ "Succinate",
            str_detect(Sample, "LMWPE") ~ "LMWPE",
            TRUE ~ NA_character_
        ),
        replicate = str_extract(Sample, "[1-3]$")
    ) %>%
    group_by(id) %>%
    t_test(
        Intensity ~ condition,
        p.adjust.method = "fdr",
        paired = TRUE
    )


volcano <- diff %>%
    left_join(ttest, by = "id") %>%
    left_join(KEGG, by = "id") %>%
    mutate(significant = case_when(
        p <= 0.05 & diff >= 1 ~ "Upregulated",
        p <= 0.05 & diff <= -1 ~ "Downregulated",
        TRUE ~ "Not significant"
    )) %>%
    mutate(protein = case_when(
        pathway == "alkane" & significant != "Not significant" ~ protein,
        pathway == "ketone" & significant != "Not significant" ~ protein,
        pathway == "plastic" & significant != "Not significant" ~ protein,
        pathway == "fatty acid" ~ NA_character_,
        TRUE ~ NA_character_
    )) %>%
    mutate(log10p = -log10(p)) %>%
    select(-c(2:7))


volcano_plot <- volcano %>%
    ggplot(aes(x = diff, y = log10p)) +
    geom_point(size = 1, aes(color = significant)) +
    geom_vline(xintercept = c(-1, 1), col = "darkgray", linetype = "dashed") +
    geom_hline(yintercept = 1.3, col = "darkgray", linetype = "dashed") +
    scale_color_manual(values = c("Upregulated" = "#ff7f00", "Downregulated" = "#7D3200", "Not significant" = "grey")) +
    geom_label_repel(aes(label = protein, fontface = "bold", segment.color = "black"), box.padding = 1, max.overlaps = Inf, min.segment.length = 0.1, size = 3.5, na.rm = TRUE) +
    labs(x = "log2(LMWPE/Succinate)", y = "-log10(p-value)") +
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
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "bottom"
    )
