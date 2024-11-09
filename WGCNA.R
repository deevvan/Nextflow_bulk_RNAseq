library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(gprofiler2)
library(GO.db)
library(igraph)
library(cowplot)
library(RColorBrewer)
library(RSkittleBrewer)
library(stringr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(gridExtra)
library(grid)
library(ggpubr)
library(WGCNA)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)


library(stats)
library(genefilter)
library(tibble)
library(ggrepel)
library(pathview)
library(AnnotationDbi)
library(biomaRt)
library(scales)
library(enrichR)
library(GOSim)
library(foreach)
library(ggridges)
library(forcats)



setwd("/path/to/working/dir")
getwd()



###################################################
# Custom function to split x or y axis tick label into two lines if longer than 30 characters
custom_labeller <- function(labels) {
  labels <- lapply(labels, function(x) {
    if (nchar(x) > 40) {
      paste0(strwrap(x, width = 55), collapse = "\n")
    } else {
      x
    }
  })
  return(labels)
}


###################################################
# function to get GO terms for clustered genes:
get_main_GO_categories <- function(gene_ids) {
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(gene = gene_ids, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  
  # Extract and return the top GO terms
  top_GO_terms <- go_enrichment[1:10, c("ID", "Description")]
  return(top_GO_terms)
}

get_main_KEGG_categories <- function(gene_ids) {
  # Perform KEGG pathway enrichment analysis
  kegg_enrichment <- enrichKEGG(gene = gene_ids, organism = "hsa")
  
  # Extract and return the top KEGG pathways
  top_KEGG_pathways <- kegg_enrichment[1:10, c("ID", "Description")]
  return(top_KEGG_pathways)
}




######################################################################################################################
#################### Enrichment plots for top up and downregulated genes for each VOC  ###############################
######################################################################################################################
######################################################################################################################


# Fetch Reactome gene sets
hs_gene2name_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::mutate(gs_name = gsub("_", " ", gs_name)) %>%
  dplyr::mutate(gs_name = str_remove(gs_name, "^REACTOME ")) %>%  
  dplyr::mutate(gs_name = tolower(gs_name))

# Fetch Reactome gene sets
hs_gene2name_GO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::mutate(gs_name = gsub("_", " ", gs_name)) %>%
  dplyr::mutate(gs_name = str_remove(gs_name, "^GOBP ")) %>%  
  dplyr::mutate(gs_name = tolower(gs_name))


# Fetch Reactome gene sets
hs_immune_genes <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") 




##############################################################################
##############################    WGCNA  #######################3#############
##############################################################################

allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Fetch Data ------------------------------------------------
data <- count_matrix1 %>%
  filter(rownames(.) %in% hs_immune_genes$ensembl_gene)

#data <- count_matrix1 

wuhan_rows <- rownames(design_matrix)[design_matrix$variant == "Wuhan"]
zero_dpi_rows <- rownames(design_matrix)[design_matrix$Timepoint == 0]

data <- data %>%
        dplyr::select(-any_of(zero_dpi_rows))



# 2. WGCNA's QC step for outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK # allOK = TRUE means all the SRR runs passed and there are no outliers

table(gsg$goodGenes) # check outliers in rows (genes)
table(gsg$goodSamples) # check outliers in columns. (SRR runs)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,] # Selecting only genes that are not outliers
gsg <- goodSamplesGenes(t(data))
table(gsg$goodGenes) # check outliers in rows (genes)

dev.off()
# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples
samples.to.be.excluded <- c('GSM7050973','GSM7050988')

data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- pheno_data_host %>% 
  filter(!Run %in% zero_dpi_rows) %>%
  filter(!Run %in% samples.to.be.excluded) %>%
  dplyr::select(-batch) %>%
  column_to_rownames("Run")


# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

#data.subset[] <- sapply(data.subset, as.numeric)
#data.subset <- round(data.subset)  # Round to nearest integer

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model

## remove all genes with counts < 15 in more than 75% of samples (77*0.75=57)
#dds75 <- dds[rowSums(counts(dds) >= 15) >= 44,] # 75% of 59
dds75 <- dds[rowSums(counts(dds) >= 15) >= 36,] # 75% of 48

nrow(dds75) # 11625 genes with counts less than 15 in over 75% of samples
# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>%
                      t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power

# Determining Maximum R square
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'signed R^2') +
  theme_classic()

# Determining Minimum mean connectivity 
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 1)

# Determined soft threshold 
soft_power <- 9 # at 9, Rsquare value was just above 80% line and mean connectivity was lowest consecutively

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

# Assigining WGCNA correlation to cor variable to avoid namespace error
temp_cor <- cor
cor <- WGCNA::cor

# Blockwise module functions's memory estimate with respect to blocksize (4GB:8k, 16GB:20k, 32GB:30k)
bwnet <- blockwiseModules(norm.counts,
                          TOMType = "signed",
                          power = soft_power,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          minModuleSize = 30,
                          mergeCutHeight = 0.25,
                          maxBlockSize = 30000,
                          deepSplit = 1,
                          verbose = 3)

cor <- temp_cor # reassign the correlation funtion back the to original cor function 


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

# get number of genes for each module
table(bwnet$colors)

# Plot dendrogram and module colors before and after merging genes
plotDendroAndColors(bwnet$dendrograms[[1]], bwnet$colors, # visualize genes before and after merging
                    c("modules"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# several colors in the unmerged modules and fewer colors in merged modules indicate that the genes belonging to 
# different unmerged modules were similar and were merged into single module after merging and merged module 
# is chosen for further analysis

# grey module = all genes that doesn't fall into other modules were assigned to the grey module





# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

# binarize categorical variables disease status
traits <- colData %>% 
  mutate(mock = ifelse(grepl('mock', variant), 1, 0)) %>% 
  dplyr::select(3)

# binarize categorical variables variants
colData$variant <- factor(colData$variant, levels = c("mock","delta","omicron"))

variant.out <- binarizeCategoricalColumns(colData$variant,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

traits <- cbind(traits, variant.out)
colnames(traits) <- c("mock", "Delta", "Omicron")


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)


# visualize module-trait association as a heatmap
# adding binary train classes to ME dataframe
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[12:14],
             y = names(heatmap.data)[1:11],
             col = c("dodgerblue2", "skyblue", "white", "pink", "red2"),
             signifSymbols = c("***", "**", "*", ""),
             signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
             corFUN = "pearson",
             main = "Module-trait plot",
             cexMain = 2,
             rotLabX = 45)

# get number of genes for each module
table(bwnet$colors)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

module.gene.mapping <- as.data.frame(bwnet$colors)

# COPY OF MAIN MODULE MAP
module.gene.mapping1 <- module.gene.mapping %>%
  tibble::rownames_to_column(var = "RowNames") %>%
  setNames(c("Gene", "Module_col"))

EG_genelist_turquoise<- module.gene.mapping %>% 
   filter(`bwnet$colors` == 'turquoise') %>% 
   rownames() %>%
   as.data.frame(.)
EG_genelist_yellow<- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames() %>%
  as.data.frame(.)

#write.csv(EG_genelist_turquoise,"turquoise_hubgenes.csv")
#write.csv(EG_genelist_yellow,"yellow_hubgenes.csv")



# Assuming module.gene.mapping1 is your dataframe with gene-module mapping
gene_lists_by_modules <- module.gene.mapping1 %>%
  group_by(Module_col) %>%
  summarize(Genes = list(unique(Gene))) %>%
  ungroup()



# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and 
# the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

# An eigengene is calculated as the first principal component (PC1) of the gene expression values within a module.
# The PC1 captures the major variation or pattern of gene expression across samples within that module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)



GS.trait_delta = as.data.frame(cor(norm.counts, traits$Delta, use = "p"))
GS.trait_omicron = as.data.frame(cor(norm.counts, traits$Omicron, use = "p"))
moduleColors <- labels2colors(bwnet$colors)
MEs = moduleEigengenes(norm.counts, moduleColors)$eigengenes
MM.trait = signedKME(norm.counts, MEs)

par(mfrow = c(1, 1))
verboseScatterplot(abs(MM.trait$kMEturquoise), abs(GS.trait_delta$V1), xlab = "Module Membership",
                     ylab = "Gene Significance for Delta", main = paste("kME.turquoise vs. GS"), col = "turquoise")
abline(h = 0.50, col = "black", lty = 2)  # Adjust color and line type as needed
abline(v = 0.80, col = "black", lty = 2)  # Adjust color and line type as needed

verboseScatterplot(abs(MM.trait$kMEyellow), abs(GS.trait_omicron$V1), xlab = "Module Membership",
                   ylab = "Gene Significance for Omicron", main = paste("kME. yellow vs. GS"), col = "yellow")
abline(h = 0.50, col = "black", lty = 2)  # Adjust color and line type as needed
abline(v = 0.80, col = "black", lty = 2)  # Adjust color and line type as needed


genes_of_interest <- cbind(GS.trait_delta,MM.trait$kMEturquoise,GS.trait_omicron,MM.trait$kMEyellow)
colnames(genes_of_interest) <- c("GS_delta","MM_turquoise","GS_omicron","MM_yellow")

# Hub genes for yellow module with delta as trait
genes_of_interest_delta <- genes_of_interest %>% 
  dplyr::select("GS_delta","MM_turquoise")
genes_of_interest_delta <- genes_of_interest_delta %>%
  filter(abs(GS_delta) >= 0.55,      
         abs(MM_turquoise) >= 0.8) 

#write.csv(genes_of_interest_delta, file = 'delta_hubs.csv', row.names = TRUE)

# Hub genes for turquoise module with omicron as trait
genes_of_interest_omicron <- genes_of_interest %>% 
  dplyr::select("GS_omicron","MM_yellow")
genes_of_interest_omicron <- genes_of_interest_omicron %>%
  filter(abs(GS_omicron) >= 0.5,      
         abs(MM_yellow) >= 0.8)    
#write.csv(genes_of_interest_omicron, file = 'omicron_hubs.csv', row.names = TRUE)

hub_genes_delta <- EG_genelist_turquoise$.
hub_genes_omicron <- EG_genelist_yellow$.




##########################################################################################################
##################################### DEG across all dpi #################################################
##########################################################################################################

design_matrix_all <- design_matrix %>% filter(variant != "Wuhan") %>% filter(Timepoint != 0)
count_matrix_all <- count_matrix1[colnames(count_matrix1) %in% rownames(design_matrix_all)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_all) %in% rownames(design_matrix_all))

# Reorder the columns of count_matrix_all to match the row names of design_matrix_all
count_matrix_all <- count_matrix_all[, rownames(design_matrix_all)]
# order check
all(colnames(count_matrix_all) == rownames(design_matrix_all))

# Converting all read counts into integers
count_matrix_all[] <- lapply(count_matrix_all, as.integer)

design_matrix_all$variant <- factor(design_matrix_all$variant, levels = c("mock","delta","omicron"))

dds_variant_all_time <- DESeqDataSetFromMatrix(countData = count_matrix_all,
                                               colData = design_matrix_all,
                                               design = ~ batch + variant)

# set the factor level
dds_variant_all_time$variant <- relevel(dds_variant_all_time$variant, ref = "mock")

deseq_variant_all_time <- DESeq(dds_variant_all_time)


# Get differential expression results
result_delta_all <- results(deseq_variant_all_time, contrast = c("variant","delta","mock"))
summary(result_delta_all)
head(result_delta_all)
result_delta_all <- as.data.frame(result_delta_all)
result_delta_all %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

result_omicron_all <- results(deseq_variant_all_time, contrast = c("variant","omicron","mock"))
summary(result_omicron_all)
head(result_omicron_all)
result_omicron_all <- as.data.frame(result_omicron_all)
result_omicron_all %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_delta_all <- result_delta_all[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta")
df_omicron_all <- result_omicron_all[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron")

log2FoldChange_table_all <- rbind(df_delta_all, df_omicron_all)

log2FoldChange_table_all <- log2FoldChange_table_all %>% filter(!is.na(log2FoldChange) & !is.na(padj)) 



log2FoldChange_table_all <- log2FoldChange_table_all %>%
  mutate(ENTREZID = mapIds(org.Hs.eg.db, log2FoldChange_table_all$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))


##########################################################################################################

### Omicron ###
diff_genes_Omicron_all <- log2FoldChange_table_all %>%
  filter(geneName %in% hub_genes_omicron) %>%
  filter(group == "Omicron") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_Delta_all <- log2FoldChange_table_all %>%
  filter(geneName %in% hub_genes_delta) %>%
  filter(group == "Delta") %>%
  filter(!is.na(ENTREZID)) 

gene_list_d_all<- diff_genes_Delta_all$log2FoldChange
names(gene_list_d_all) <- diff_genes_Delta_all$ENTREZID
gene_list_d_all= sort(gene_list_d_all, decreasing = TRUE)


gene_list_o_all<- diff_genes_Omicron_all$log2FoldChange
names(gene_list_o_all) <- diff_genes_Omicron_all$ENTREZID
gene_list_o_all= sort(gene_list_o_all, decreasing = TRUE)


# ClusterprofilR GSEA GO pathway
gse_d_all <- clusterProfiler::GSEA(gene_list_d_all,
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   minGSSize = 15,
                                   gson = NULL,
                                   TERM2GENE = hs_gene2name_GO)

gse_o_all <- clusterProfiler::GSEA(gene_list_o_all,
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   minGSSize = 15,
                                   gson = NULL,
                                   TERM2GENE = hs_gene2name_GO)


GOplot_d_all <- ggplot(data = gse_d_all@result%>%
                         top_n(15, rank) %>%
                         mutate(Description = fct_reorder(Description, setSize)),  # Reorder y-axis by NES
                       mapping = aes(x = NES, y = Description, 
                                     color = p.adjust, size = setSize)) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("GO for Delta Hub genes") +
  scale_y_discrete(labels = custom_labeller)
GOplot_d_all



GOplot_o_all <- ggplot(data = gse_o_all@result %>%
                         top_n(15, rank) %>%
                         mutate(Description = fct_reorder(Description, setSize)),  # Reorder y-axis by NES
                       mapping = aes(x = NES, y = Description, 
                                     color = p.adjust, size = setSize)) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("GO for Omicron Hub genes") +
  scale_y_discrete(labels = custom_labeller)
GOplot_o_all



















