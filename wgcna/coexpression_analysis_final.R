####Script to perform Wighted Gene Co-Expression Analysis (WGCNA) and plot Figure 6

# Author: Dave Edison Rojas Calderon

# # Install / Load WGCNA
# install.packages("WGCNA")
library(WGCNA)
# Enable multiprocessing
enableWGCNAThreads()
# # Load tidyverse
# install.packages("tidyverse")
library(tidyverse)
# # Load flashClust
# install.packages("flashClust")
library(flashClust)
# Load ggplot2
library(ggplot2)
# # Load ggpubr
# install.packages("ggpubr")
library(ggpubr)
# Load openxlsx
library(openxlsx)
# Load upsetR
library(UpSetR)
# Load gridG
library(grid)
# Load cowplot
library(cowplot)

# DefinegridGraphics# Define functions
# ------------------------------------------------------------------------------
# Function to check if a directory exist and create it when it does not
check_create_dir = function(directory){
  # Check if the given directory exists
  if (file.exists(directory)) {
    print(paste0("The directory already exists: ", directory))
  } else{
    # Create directory when it does not exist
    dir.create(directory)
    print(paste0("The directory has been created: ", directory))
  }
}

# Function to build custom gene sets (using GeneID instead of KO)
get_kegg_pathway_sets = function(df, verbose = TRUE) {
  # Check expected columns
  if (!all(c("GeneID", "ko_id") %in% names(df))) {
    stop("Input data frame must have columns: 'GeneID' and 'ko_id'")
  }
  
  # Remove rows with missing KO
  df = df %>% filter(!is.na(ko_id), ko_id != "")
  
  # Make sure we only query unique KOs
  unique_kos = unique(df$ko_id)
  
  custom_gene_sets = list()
  
  for (ko in unique_kos) {
    if (verbose) message("Processing ", ko)
    
    # Safely query KEGG
    res = tryCatch({
      keggGet(paste0("ko:", ko))
    }, error = function(e) NULL)
    
    # Skip if no result
    if (is.null(res) || length(res) == 0) next
    
    pathways = res[[1]]$PATHWAY
    if (is.null(pathways)) next
    
    # Extract pathway names (values of the named vector)
    pathway_names = unname(pathways)
    
    # Find all GeneIDs associated with this KO
    genes_for_ko = df %>%
      filter(ko_id == ko) %>%
      pull(GeneID) %>%
      unique()
    
    # Add to each pathway in the list
    for (pw in pathway_names) {
      if (pw %in% names(custom_gene_sets)) {
        custom_gene_sets[[pw]] = unique(c(custom_gene_sets[[pw]], genes_for_ko))
      } else {
        custom_gene_sets[[pw]] = genes_for_ko
      }
    }
  }
  
  return(custom_gene_sets)
}

# Loading input data
# ------------------------------------------------------------------------------
# Defining input paths
proteomics_path = file.path("data/RhRMS-31_intensity_imputed.tsv")
dea_res_path = file.path("data/Table S1_V1.xlsx")
dram_path = file.path("data/annotations.tsv")

# Define output directories
outdir_plot = file.path("plots/")
outdir_tables = file.path("data/proc_data/")
# Create directories if they do not exist
check_create_dir(outdir_plot)
check_create_dir(outdir_tables)

# Load data
proteomics = read_delim(proteomics_path,
                        delim = "\t") %>% 
  # Clean Protein names
  mutate(Protein = str_replace(Protein, "2_01_Rhodococcus_sp002259485", "RhASF")) %>% 
  # Rename Protein names
  rename(GeneID = Protein)

# Load Annotation DRAM
ko_annotation = read_delim(dram_path,
                           delim = "\t")  %>% 
  # Rename first column
  rename(GeneID = `...1`) %>% 
  # Correct Isolate naming
  mutate(GeneID = str_replace(GeneID, "2_01_Rhodococcus_sp002259485", "RhASF")) %>% 
  # Keep interesting columns
  select(GeneID, ko_id)

# Load DEA results
dea_res = readxl::read_excel(dea_res_path,
                             sheet = "Table S1A",
                             range = "B4:N3620") %>% 
  # Select only interesting columns
  select(`Gene ID`, Protein, Pathway) %>% 
  # Rename columns with problematic names
  rename(GeneID = `Gene ID`) %>% 
  mutate(GeneID = str_replace(GeneID, "RhRMS", "RhASF")) 

# Save annotation
annotation = proteomics %>% 
  select(GeneID, Annotation) %>% 
  # Add ko_id
  left_join(ko_annotation,
            by = "GeneID") %>% 
  left_join(dea_res,
            by = "GeneID")
  
# Data as matrix
expr = proteomics %>% 
  select(-Annotation) %>% 
  column_to_rownames(var = "GeneID") %>% 
  as.matrix()

# Define Substrate_df
Substrate_df = data.frame(Sample = colnames(expr)) %>% 
  mutate(Substrate = ifelse(str_detect(Sample, "Succinate"), "Succinate", "PELW"))

# Looking for outliers
# ------------------------------------------------------------------------------
# Transposet expression data
expr_t = expr %>% 
  t()

# Make a hierachical cluerting of sample to check for outlieras
sample_tree = hclust(dist(expr_t), method = "average")

# Define plot out path
plot_path = file.path(outdir_plot, "sample_clustering_outliers.pdf")
# Initialize pdf
pdf(plot_path)

# Plot the dendogram
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sample_tree, main = "Sample clustering to detect outliers",
     sub="",
     xlab="Samples",
     cex.lab = 1.5, 
     cex.axis = 1.5,
     cex.main = 2)
# Add the cut-off
abline(h = 80, col = "red")

# Samples are clustered in two groups, succinate and PELW, as expected (no outliers).

# Prepare annotation data for posterior cluster association analysis
# ------------------------------------------------------------------------------
# Get sample names
samples = colnames(expr)

# Define colors for Substrates
Substrate_colors = c("grey", "grey", "grey", "#8B2323", "#8B2323", "#8B2323")

# Plot a sample dendogram with the colors below
plotDendroAndColors(sample_tree, Substrate_colors,
                    groupLabels = "Substrate",
                    main = "Sample dendrogram and trait heatmap")

# close pdf
dev.off()

# Creating Network
# ------------------------------------------------------------------------------
# Choose a set of sft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  expr_t,
  powerVector = powers,
  verbose = 5
)

# Get results as df
res = sft$fitIndices %>% 
  as.data.frame()

# Index the scale free topology adjust as a function of the power soft thresholding
sft_plot = res %>% 
  ggplot(aes(x = Power,
           y = -sign(slope) * SFT.R.sq,
           label = Power)) +
  geom_label() +
  theme_bw() +
  labs(x = "Soft Threshold (power)",
       y = "Scale Free Topology Model Fit,signed R^2") +
  geom_hline(yintercept = 0.9,
             color = "red")

# Show
plot(sft_plot)
# Save
ggsave(sft_plot,
       path = outdir_plot,
       filename = "soft_threshold.pdf",
       device = "pdf")

# Connectivity mean as a function of soft power thresholding
mean_k_plot = res %>% 
  ggplot(aes(x = Power,
             y =  mean.k.,
             label = Power)) +
  geom_label() +
  theme_bw() +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity")

# Show
plot(mean_k_plot)
# Save
ggsave(mean_k_plot,
       path = outdir_plot,
       filename = "mean_connectivity.pdf",
       device = "pdf")

# Define picked power
picked_power = 20
# Calculate adjacency matrix
adjacency = adjacency(expr_t,
                      power = picked_power,
                      type = "signed")
# Transform adjacency matrix in a topological overlap
tom = TOMsimilarity(adjacency)
# Add names to tom
colnames(tom) = colnames(adjacency)
rownames(tom) = rownames(adjacency)
# Calculate dissimilarity
diss_tom = 1 - tom

# Hierarchical Clustering Analysis to determine modules
# ------------------------------------------------------------------------------
# Create the dendrogram 
gene_tree = hclust(as.dist(diss_tom), method = "average")
# Determine modules
modules = cutreeDynamic(dendro = gene_tree,
                        distM = diss_tom,
                        deepSplit = 2,
                        pamRespectsDendro = FALSE,
                        minClusterSize = 30)
table(modules) # Observe modules

# Assign colors to modules to add it below the dendogram
module_colors = labels2colors(modules)
table(module_colors) # Visualize colors of modules

# # Count total number of proteins for each module
# module_summary = tibble(color = module_colors) %>%
#   count(color, name = "occurrences") %>%
#   mutate(module = paste0("ME", color)) %>%
#   select(module, occurrences) %>% 
#   mutate(module = factor(module, levels = module))
# # Plot number of proteins per module
# prots_module_plot = module_summary %>% 
#   ggplot(aes(x = module,
#              y = occurrences,
#              fill = module)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = gsub("ME", "", module_summary$module)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Module",
#        y = "Total number",
#        title = "Total number of Proteins per Module")
# 
# # Show
# plot(prots_module_plot)
# # Save
# ggsave(prots_module_plot,
#        path = outdir_plot,
#        filename = "protein_per_module.pdf",
#        device = "pdf")
# 
# # Define plot path
# plot_path = file.path(outdir_plot, "protein_dendogram_modules.pdf")
# # Initialize pdf
# pdf(plot_path)
# # Plots the protein dendrogram with the module colors
# plotDendroAndColors(gene_tree, module_colors, "Module",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Protein dendrogram and module colors")
# # Close pdf
# dev.off()

# Grouping (Clustering) genes into modules
# ------------------------------------------------------------------------------
# Calculate eigengenes
me_list = moduleEigengenes(expr = expr_t,
                           colors = module_colors)
mes = me_list$eigengenes

# Calculating module dissimilarity
me_diss = 1 - cor(mes)

# Clustering the eigengenes modules
me_tree = hclust(as.dist(me_diss), method = "average")

# Define plot path
plot_path = file.path(outdir_plot, "dendogram_modules_to_merge.pdf")
# Initialize pdf
pdf(plot_path)
# Set mergin sizes
par(mar = c(0,4,2,0))
# Scale the graphic
par(cex = 0.6);
# Plot dendogram
plot(me_tree)
# Add line at height 0.25 (correlation of 0.75)
abline(h=.25, col = "red")
# Close pdf
dev.off()

# Merge modules
merge = mergeCloseModules(expr_t, module_colors, cutHeight = .25)
# The merged module colors, assigning one color to each module
merged_colors = merge$colors
# Eigengenes of the new merged modules
merged_mes = merge$newMEs

# Count total number of proteins for each module
module_merged_summary = tibble(color = merged_colors) %>%
  count(color, name = "occurrences") %>%
  mutate(module = color) %>%
  select(module, occurrences)
# Get order to plot Proteins
merged_module_order = module_merged_summary %>% 
  arrange(desc(occurrences)) %>% 
  pull(module)

# Plot number of proteins per module
prots_merged_module_plot = module_merged_summary %>% 
  # Add correct level order to modules
  mutate(module = factor(module, levels = merged_module_order)) %>% 
  ggplot(aes(x = module,
             y = occurrences,
             fill = module)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = merged_module_order) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "",
       y = "Total number of Proteins") +
  theme(legend.position = "none")

# Module and Traits matching
# ------------------------------------------------------------------------------
# Calculate the association between expression profile and trait by using the
# eigengene gene significance (correlation between trait and module eigengene).
# Define numbers of genes and samples
n_genes = ncol(expr_t)
n_samples = nrow(expr_t)

# Get eigengenes per cluster
mes_0 = moduleEigengenes(expr_t, merged_colors)$eigengenes
# Reorder modules so similar modules are next to each other
mes_0 = orderMEs(mes_0)
module_order = names(mes_0) %>%
  gsub("ME","", .)

# Add Substrate names
mes_0$Substrate = row.names(mes_0)

# tidy & plot data
m_me = mes_0 %>%
  pivot_longer(-Substrate) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

# Pepare data for Complex Heatmap
me_mat = m_me %>% 
  pivot_wider(names_from = "Substrate",
              values_from = "value") %>% 
  column_to_rownames(var = "name") %>% 
  as.matrix()

# Module trait-relationship (at trait level not sample)
# ------------------------------------------------------------------------------
# Convert your Substrate_df to numeric dummy variable
Substrate_numeric = Substrate_df %>%
  column_to_rownames("Sample") %>%
  rename(Substrate = Substrate) %>% 
  mutate(Substrate = ifelse(Substrate == "Succinate", 0, 1))

# Calculate correlation
module_trait_corr = cor(t(me_mat), Substrate_numeric, use = "p")
# Define a list to store the results
me_pval = c()

# Initialize a loop to run a t.test between Succinate and PELW
for (i in 1:ncol(merge$newMEs)) {
  me_pval[i] = t.test(merge$newMEs[str_detect(rownames(merge$newMEs), "Succinate"), i], merge$newMEs[!str_detect(rownames(merge$newMEs), "Succinate"), i])$p.value
}
# Table with module names, p-values, adjusted p-values, scaled "significance"
# by computing -log10(adjusted p-values)
me_sig_tbl = tibble(Modules = str_replace(colnames(merge$newMEs), "ME", ""),
                    MEpval = me_pval,
                    MEpval.adj = p.adjust(me_pval, method = "fdr"),
                    MEsignificance = -log10(MEpval.adj)) %>% 
  arrange(Modules)
# Get module-trait adjusted p-value
module_trait_pval = data.frame(Substrate_padj = me_sig_tbl$MEpval.adj,
                               Module = me_sig_tbl$Modules) %>%
  column_to_rownames(var = "Module") %>% 
  as.matrix()
# Order rows
module_trait_pval = module_trait_pval[rownames(module_trait_corr), , drop = FALSE]

# Plot
# Define plot path
plot_path = file.path(outdir_plot, "module_trait_relationship.pdf")
# Initialize pdf
pdf(plot_path,
    width = 6)
# Will display correlations and their p-values
textMatrix = paste(signif(module_trait_corr, 2), "\n(",
                   signif(module_trait_pval, 1), ")", sep = "");
dim(textMatrix) = dim(module_trait_corr)
par(mar = c(6, 8.5, 3, 1))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module_trait_corr,
               xLabels = names(Substrate_numeric),
               yLabels = names(merged_mes),
               ySymbols = rownames(module_trait_corr),
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationship"),
               yColorLabels = TRUE)
# Close pdf
dev.off()

# Module Eigengene Significance
# ------------------------------------------------------------------------------
# Plot significance for each module
me_sig_plot = me_sig_tbl %>% 
  ggplot(aes(x = Modules, y = MEsignificance, fill = Modules)) +
  geom_col(col="black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  scale_fill_manual(values=gsub("ME", "", module_merged_summary$module)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
        legend.position="none") +
  ylab("Module eigengene significance")

# Show
plot(me_sig_plot)
# Save
ggsave(me_sig_plot,
       path = outdir_plot,
       filename = "ME_significance.pdf",
       device = "pdf")

# Module green and blue are the only significant module, being green the most interesting for our study

# Make a df with module information
annotation_updated = data.frame(GeneID = annotation$GeneID,
                            Annotation = annotation$Annotation,
                            Protein = annotation$Protein,
                            Pathway = annotation$Pathway,
                            Module = merge$colors)

# Recover pathways for proteins with KO ids
# ------------------------------------------------------------------------------
# Build the custom gene sets (Takes a while)
# custom_gene_sets = get_kegg_pathway_sets(annotation)
# save(custom_gene_sets, file = file.path("data/custom_gene_sets.RDa"))
load(file.path("data/custom_gene_sets.RDa"))

# Intersection analysis of DE and co-expressed proteins in module turquoise
# ------------------------------------------------------------------------------
# Define proteins involved in biofilm formation
biofilm_prots = c(custom_gene_sets[["Biofilm formation - Vibrio cholerae"]],
                  custom_gene_sets[["Biofilm formation - Escherichia coli"]],
                  custom_gene_sets[["Biofilm formation - Pseudomonas aeruginosa"]])
# Define proteins in Fatty Acid Metabolism
fam_prots = c(custom_gene_sets[["Fatty acid metabolism"]])
# Define descriptions involved in biofilm formation
biofilm_descr = c("exopolysaccharide production protein ExoY (db=kegg)",
                  "GDPmannose 4,6-dehydratase [EC:4.2.1.47] (db=kegg)",
                  "GDP-L-fucose synthase [EC:1.1.1.271] (db=kegg)",
                  "putative colanic acid biosynthesis glycosyltransferase WcaI (db=kegg)",
                  "Glycosyltransferase Family 4 [PF13439.11]; Glycosyl transferases group 1 [PF00534.25]; Glycosyl transferase 4-like domain [PF13579.11] (db=pfam)",
                  "Polysaccharide biosynthesis protein [PF13440.11] (db=pfam)")
# Define descriptions for surfactant production
surfactant_descr = c(
  "diacylglycerol O-acyltransferase / trehalose O-mycolyltransferase / mycolyltransferase Ag85 [EC:2.3.1.20 2.3.1.122 2.3.1.-] (db=kegg)",
  "diacylglycerol O-acyltransferase / trehalose O-mycolyltransferase / mycolyltransferase Ag85 [EC:2.3.1.20 2.3.1.122 2.3.1.-] (db=kegg)"
)
# Add the new interesting pathways to the annotation
annotation_updated = annotation_updated %>% 
  # Add Pathways
  mutate(Annotation = trimws(Annotation),
         Pathway = ifelse(GeneID %in% biofilm_prots, "Biofilm formation", Pathway),
         Pathway = ifelse(Annotation %in% biofilm_descr, "Biofilm formation", Pathway),
         Pathway = ifelse(GeneID %in% fam_prots & is.na(Pathway), "Fatty Acid Metabolism", Pathway),
         Pathway = ifelse(Annotation == "diacylglycerol O-acyltransferase / trehalose O-mycolyltransferase / mycolyltransferase Ag85 [EC:2.3.1.20 2.3.1.122 2.3.1.-] (db=kegg)", "Surfactant production", Pathway))
# Extract the unique pathways for upset plot
biofilm = annotation_updated %>% 
  filter(Pathway == "Biofilm formation") %>% 
  pull(GeneID)
fats = annotation_updated %>% 
  filter(Pathway == "fatty acid") %>% 
  pull(GeneID)
fam = annotation_updated %>% 
  filter(Pathway == "Fatty Acid Metabolism") %>% 
  pull(GeneID)
alkane = annotation_updated %>% 
  filter(Pathway == "alkane") %>% 
  pull(GeneID)
surfactant = annotation_updated %>% 
  filter(Pathway == "Surfactant production") %>% 
  pull(GeneID)
ketone = annotation_updated %>% 
  filter(Pathway == "ketone") %>% 
  pull(GeneID)
plastic = annotation_updated %>% 
  filter(Pathway == "plastic") %>% 
  pull(GeneID)

# Extract Genes in Module Green
green = annotation_updated %>% 
  filter(Module == "green") %>% 
  pull(GeneID)
# Define gene sets to intersect
gene_sets = list(
  Green     = green,
  Biofilm_Formation = biofilm,
  FATS = fats,
  FAM = fam,
  Alkane_Degradation = alkane,
  Surfactant_Production = surfactant,
  Ketone_Degradation = ketone,
  Plastic_Degradation = plastic
)
# Give format required by UpSetR
upset_input = fromList(gene_sets)
# Assign gene names as rownames
rownames(upset_input) = unique(unlist(gene_sets))

# Make Final Figure
# ------------------------------------------------------------------------------
# Save as pdf
png(file.path("plots/panel_A.png"))
# Panel A: Cluster Dendogram
plotDendroAndColors(
      gene_tree,
      merged_colors,
      "Modules",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05,
      main = NA
    )
dev.off()

# Panel B: Number of Proteins per Module
ggsave(plot = prots_merged_module_plot,
       filename = "panel_B.png",
       path = file.path("plots/"),
       device = "png",
       height = 4,
       width = 4)

# Save as pdf
png(file.path("plots/panel_C.png"))
# Define colors
orangeWhiteBrown = colorRampPalette(c("#ff7f00", "#FFFFFF", "#7D3200"))
# Display correlations and their p-values
textMatrix_me = paste(signif(me_mat, 2), "\n(",
                      signif(module_trait_pval, 1), ")", sep = "")
# Panel C: Module Eigengene per Sample
labeledHeatmap(Matrix = me_mat,
               xLabels = colnames(me_mat),
               yLabels = paste0("ME",rownames(me_mat)),
               ySymbols = rownames(me_mat),
               colors = rev(orangeWhiteBrown(50)),
               textMatrix = textMatrix_me,
               setStdMargins = TRUE,
               cex.text = 1,
               main = NA,
               yColorLabels = TRUE,
               verticalSeparator.x = seq(1, 5, 1),
               horizontalSeparator.y = seq(1, 9, 1),
               verticalSeparator.col = "#FFFFFF",
               horizontalSeparator.col = "#FFFFFF",
               verticalSeparator.lwd = 2,
               horizontalSeparator.lwd = 2,
               verticalSeparator.ext = 0,
               horizontalSeparator.ext = 0)
dev.off()


# Panel D: UpSet Plot of Up-regulated genes


upset(
  upset_input,
  sets = colnames(upset_input),
  intersections = list(
    list("Green", "Alkane_Degradation"),
    list("Green", "Biofilm_Formation"),
    list("Green", "FATS"),
    list("Green", "FAM"),
    list("Green", "Ketone_Degradation"),
    list("Green", "Plastic_Degradation"),
    list("Green", "Surfactant_Production")
  ),
  order.by = "freq",    # order by intersection size
  keep.order = FALSE,     # keep set order as defined
  queries = list(
    list(query = intersects,
         params = list("Alkane_Degradation", "Green"),
         color = "#ff7f00", active = TRUE),
    list(query = intersects,
         params = list("Biofilm_Formation", "Green"),
         color = "#ff7f00", active = TRUE),
    list(query = intersects,
         params = list("FATS", "Green"),
         color = "#ff7f00", active = TRUE),
    list(query = intersects,
         params = list("FAM", "Green"),
         color = "#ff7f00", active = TRUE),
    list(query = intersects,
         params = list("Ketone_Degradation", "Green"),
         color = "#ff7f00", active = TRUE),
    list(query = intersects,
         params = list("Plastic_Degradation", "Green"),
         color = "#ff7f00", active = TRUE),
    list(query = intersects,
         params = list("Surfactant_Production", "Green"),
         color = "#ff7f00", active = TRUE)
  ),
  shade.color = "#969696",
  sets.bar.color = "#7D3200",
  set_size.show = TRUE,
  text.scale = c(1.5, 1.5, 1.5, 1, 1.5, 2),
  # 👇 Add this line
  set_size.scale_max = 1900
)
panel_d = recordPlot()

# Save as pdf
png(file.path("plots/panel_D.png"))
replayPlot(panel_d)
dev.off()

# Paths to your four PDF panels
files = c("plots/panel_A.png",
          "plots/panel_C.png",
          "plots/panel_B.png",
          "plots/panel_D.png")

render_pdf_to_grob = function(path) {
  # Read back as grob
  img = png::readPNG(path)
  panel = rasterGrob(img, interpolate = TRUE)
}


grobs = lapply(files, render_pdf_to_grob)

# Combine 2×2 grid
multi_panel = plot_grid(
  grobs[[1]], grobs[[2]],
  grobs[[3]], grobs[[4]],
  align = "hv",
  labels = c("A", "C", "B", "D"),
  scale = c(1, 1, 1, 1),
  ncol = 2,
  label_size = 20
)

# Save final figure
ggsave("plots/combined_panels.pdf", multi_panel, width = 12, height = 10)

# Save Tables
# ------------------------------------------------------------------------------
# Create workbook
wb = createWorkbook()

# Add sheet and write Proteins and their MNodules
addWorksheet(wb, "Protein Modules")
writeData(wb, "Protein Modules", annotation_updated)

# Add sheet and write Module eigengenes
addWorksheet(wb, "Module Eigengenes")
writeData(wb, "Module Eigengenes", me_mat, rowNames = TRUE)

# Add sheet and write Module significance
addWorksheet(wb, "Module-Trait Significance")
writeData(wb, "Module-Trait Significance", me_sig_tbl)

# Add and write Module-Trait correlation
addWorksheet(wb, "Module-Trait Correlation")
writeData(wb, "Module-Trait Correlation", module_trait_corr, rowNames = TRUE)

# Add and write number of proteins per module
addWorksheet(wb, "Proteins per module")
writeData(wb, "Proteins per module", module_merged_summary)

# Add and write Gene Set intersections
addWorksheet(wb, "Gene Set Intersections")
writeData(wb, "Gene Set Intersections", upset_input, rowNames = TRUE)

# Define path to save
excel_path = file.path("data/proc_data/Table_S3.xlsx")
# Save the workbook
saveWorkbook(wb, excel_path, overwrite = TRUE)
