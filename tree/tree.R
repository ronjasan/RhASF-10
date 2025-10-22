#### Script to make phylogenetic tree plot Figure 1

# Author: Ronja Marlonsdotter Sandholm and Dave Edison Rojas Calderon

library(tidyverse)
library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(ape)
library(ggplot2)
library(aplot)
library(ggstar)
library(cowplot)
library(readxl)


unrooted_tree <- read.tree("data/gtdbtk.unrooted.tree")
rooted_tree <- root(unrooted_tree,
    outgroup = "Corynebacterium_antarcticum_CCM_8835_GCA_016595285.2",
    resolve.root = TRUE
)

checkm_res <- read_csv("data/checkm_output.csv") %>%
    select(Bin_Id, Completeness, Contamination) %>%
    rename(ID = Bin_Id)

info <- read_xlsx("data/strain_info.xlsx") %>%
    right_join(checkm_res) %>%
    mutate(Source = fct_relevel(Source, "Host: animal", "Host: plant", "Soil", "Water", "Other")) %>%
    mutate(Plastic = fct_relevel(Plastic, "Polyethylene", "Crude oil", "Alkane", "No")) %>%
    mutate(Genome_size = Genome_size / 1e6)

d_dummy <- data.frame(
    label = rooted_tree$tip.label,
    category = as.character(1:length(rooted_tree$tip.label))
)
p_dummy <- ggplot(data = d_dummy, aes(y = label, x = category)) +
    theme_void()

info_long <- info %>%
    select(ID, Completeness, Contamination) %>%
    mutate(Completeness = Completeness - Contamination) %>%
    dplyr::rename(`Completeness (%)` = Completeness, `Contamination (%)` = Contamination) %>%
    pivot_longer(cols = c(`Completeness (%)`, `Contamination (%)`), names_to = "Metric", values_to = "Value")

ggt <- ggtree(rooted_tree, aes(color = Plastic)) %<+% info %<+% info_long +
    geom_tiplab(aes(label = Species),
        align = TRUE,
        offset = 0.001,
        size = 5
    ) +
    geom_star(aes(starshape = Source, fill = Source), size = 3, color = "black") +
    scale_shape_manual(values = c(28, 23, 11, 13, 15, 1)) +
    scale_fill_manual(values = c("#7D3200", "#B78560", "#ff7f00", "#fdae6b", "#d9d9d9")) +
    scale_color_manual(values = c("#ff7f00", "#7D3200", "#969696", "#231F20")) +
    # geom_treescale()+
    xlim(0, 0.7) +
    theme_tree2() +
    theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        axis.line = element_blank()
    )

mimag <- ggplot(filter(info_long), aes(x = Value, y = ID, fill = Metric)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    scale_fill_manual(values = c(microshades_palettes$micro_orange[5], "#7D3200")) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black", linewidth = 0.75),
        axis.line.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 12, ),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)
    ) +
    labs(
        fill = "Genome quality"
    )

size <- ggplot(filter(info), aes(x = Genome_size, y = ID)) +
    geom_bar(stat = "identity", position = "stack", fill = "#969696", width = 0.75) +
    scale_fill_manual(values = microshades_palettes$micro_gray[3]) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black", linewidth = 0.75),
        axis.line.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 12, angle = 90, vjust = 0.5, hjust = 1)
    ) +
    labs(
        fill = "Genome size (Mb)"
    )

tree <- p_dummy %>%
    insert_right(ggt + theme(legend.position = "none"), width = 30) %>%
    insert_right(mimag + theme(legend.position = "none"), width = 3) %>%
    insert_right(size + theme(legend.position = "none"), width = 3)
