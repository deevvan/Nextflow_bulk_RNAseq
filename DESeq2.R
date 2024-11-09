# BiocManager::install("ComplexHeatmap")

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



setwd("/mmfs1/projects/changhui.yan/DeewanB/VOC_SRR/deseq2_human/backup_wgcna")
getwd()

# Load the saved R session
load(file = "deseq2_human_rsemcount.RData")



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





##############################################################################################################################
##############################################################################################################################
################################################## HUMAN WITH HEALTHY CONTROLS ###############################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################



metadata_ncbi <- read.csv("SRP427025_metadata_from_ncbi.csv")
pheno_data_host <- metadata_ncbi[,c("Sample.Name","virus","batch","Timepoint")]
pheno_data_host <- pheno_data_host %>%
  filter(virus != "SA")
colnames(pheno_data_host) <- c("Run","variant","batch","Timepoint")


# Load the read count matrix
count_matrix <- read.csv("SRP427025_rsem_count_matrix.csv", header=TRUE, row.names=1)


# ## Converting Ensembl IDs (eg: ENSG00000000003) into gene IDs
# library(biomaRt)
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# ensembl_ids <- rownames(count_matrix)
# 
# # Query Ensembl for gene symbols based on Ensembl gene IDs
# gene_info <- getBM(
#   attributes = c("ensembl_gene_id", "hgnc_symbol"),  # hgnc_symbol gives the gene symbol
#   filters = "ensembl_gene_id",
#   values = ensembl_ids,
#   mart = ensembl
# )
# 
# 
# # Merge gene symbols into the count matrix
# count_matrix$ensembl_gene_id <- rownames(count_matrix)
# # Perform the merge
# count_matrix <- merge(
#   count_matrix,
#   gene_info,
#   by.x = "ensembl_gene_id",
#   by.y = "ensembl_gene_id",
#   all.x = TRUE
# )
# # Filter out rows where hgnc_symbol is blank or NA
# count_matrix <- count_matrix %>%
#   filter(hgnc_symbol != "" & !is.na(hgnc_symbol))%>%
#   distinct(hgnc_symbol, .keep_all = TRUE)
# 
# # Making sure all gene names are unique
# table(duplicated(count_matrix$hgnc_symbol))
# 
# # Convert the 'hgnc_symbol' column to row names
# count_matrix1 <- count_matrix %>%
#   column_to_rownames(var = "hgnc_symbol") %>%
#   dplyr::select(-ensembl_gene_id)
# 
# count_matrix1 <- as.data.frame(count_matrix1)

count_matrix1 <- as.data.frame(count_matrix)


# Subset pheno_data_host based on the matched column names
pheno_data_host <- pheno_data_host %>% filter(variant != "Wuhan" & Timepoint != 0)

pheno_data_host <- pheno_data_host %>% filter(pheno_data_host$Run %in% 
                                                intersect(pheno_data_host$Run, colnames(count_matrix1)))


# Subsetting count matrix to only contain the SRR runs that were selected in metadata file
count_matrix1 <- count_matrix1[, intersect(colnames(count_matrix1), pheno_data_host$Run)]


count_matrix1[is.na(count_matrix1)] <- 0

# Design Matrix 
design_matrix <- as.data.frame(cbind(pheno_data_host$Run,
                                     pheno_data_host$variant,
                                     pheno_data_host$batch,
                                     pheno_data_host$Timepoint))
design_matrix %<>% {rownames(.) <- .$V1; .[, -1]} %>% 
                    {colnames(.) <- c("variant","batch","Timepoint"); .} 

design_matrix_out <- t(design_matrix)

#################################################

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix1) %in% rownames(design_matrix))

# Reorder the columns of count_matrix to match the row names of design_matrix
count_matrix1 <- count_matrix1[, rownames(design_matrix)]
# order check
all(colnames(count_matrix1) == rownames(design_matrix))


# Converting all read counts into integers
count_matrix1[] <- lapply(count_matrix1, as.integer)

table(design_matrix$Timepoint)


design_matrix$variant <- factor(design_matrix$variant, levels = c("mock","delta","omicron"))
design_matrix$batch <- factor(design_matrix$batch)
design_matrix$Timepoint <- factor(design_matrix$Timepoint)





######################################################################################################
##################################### DEG at each hpi#################################################
######################################################################################################

######### Expression at each timepoint:######### 

########### 1 hours ############
design_matrix_1 <- filter(design_matrix, Timepoint %in% c(1))

count_matrix_1 <- count_matrix1[colnames(count_matrix1) %in% rownames(design_matrix_1)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_1) %in% rownames(design_matrix_1))

# Reorder the columns of count_matrix_1 to match the row names of design_matrix_1
count_matrix_1 <- count_matrix_1[, rownames(design_matrix_1)]
# order check
all(colnames(count_matrix_1) == rownames(design_matrix_1))

# Converting all read counts into integers
count_matrix_1[] <- lapply(count_matrix_1, as.integer)

design_matrix_1$variant <- factor(design_matrix_1$variant, levels = c("mock","delta","omicron"))
design_matrix_1$batch <- factor(design_matrix_1$batch)

dds_variant_time_1 <- DESeqDataSetFromMatrix(countData = count_matrix_1,
                                             colData = design_matrix_1,
                                             design = ~ batch + variant)

# set the factor level
dds_variant_time_1$variant <- relevel(dds_variant_time_1$variant, ref = "mock")


# DESeq on 1h
deseq_variant_time_1 <- DESeq(dds_variant_time_1)


# Get differential expression results
result_delta_1 <- results(deseq_variant_time_1, contrast = c("variant","delta","mock"))
summary(result_delta_1)
head(result_delta_1)
result_delta_1 <- as.data.frame(result_delta_1)
result_delta_1 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .} # rownames written to geneName column and rownames removed


result_omicron_1 <- results(deseq_variant_time_1, contrast = c("variant","omicron","mock"))
summary(result_omicron_1)
head(result_omicron_1)
result_omicron_1 <- as.data.frame(result_omicron_1)
result_omicron_1 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_delta_1 <- result_delta_1[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_1")
df_omicron_1 <- result_omicron_1[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_1")

log2FoldChange_table_1 <- rbind(df_delta_1, df_omicron_1)

log2FoldChange_table_1 <- log2FoldChange_table_1 %>% filter(!is.na(log2FoldChange) & !is.na(padj))


# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_1 <- log2FoldChange_table_1%>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05)


# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
diff_genes_1_immune <- diff_genes_1[diff_genes_1$geneName %in% hs_immune_genes$ensembl_gene, ]






########### 2 hours ############
design_matrix_2 <- filter(design_matrix, Timepoint %in% c(2))

count_matrix_2 <- count_matrix1[colnames(count_matrix1) %in% rownames(design_matrix_2)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_2) %in% rownames(design_matrix_2))

# Reorder the columns of count_matrix_2 to match the row names of design_matrix_2
count_matrix_2 <- count_matrix_2[, rownames(design_matrix_2)]
# order check
all(colnames(count_matrix_2) == rownames(design_matrix_2))

# Converting all read counts into integers
count_matrix_2[] <- lapply(count_matrix_2, as.integer)

design_matrix_2$variant <- factor(design_matrix_2$variant, levels = c("mock","delta","omicron"))
design_matrix_2$batch <- factor(design_matrix_2$batch)

dds_variant_time_2 <- DESeqDataSetFromMatrix(countData = count_matrix_2,
                                             colData = design_matrix_2,
                                             design = ~ batch + variant)

# set the factor level
dds_variant_time_2$variant <- relevel(dds_variant_time_2$variant, ref = "mock")
dds_variant_time_2$variant


# DESeq on 1h
deseq_variant_time_2 <- DESeq(dds_variant_time_2)


# Get differential expression results
result_delta_2 <- results(deseq_variant_time_2, contrast = c("variant","delta","mock"))
summary(result_delta_2)
head(result_delta_2)
result_delta_2 <- as.data.frame(result_delta_2)
result_delta_2 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_omicron_2 <- results(deseq_variant_time_2, contrast = c("variant","omicron","mock"))
summary(result_omicron_2)
head(result_omicron_2)
result_omicron_2 <- as.data.frame(result_omicron_2)
result_omicron_2 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_delta_2 <- result_delta_2[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_2")
df_omicron_2 <- result_omicron_2[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_2")

log2FoldChange_table_2 <- rbind(df_delta_2, df_omicron_2)

log2FoldChange_table_2 <- log2FoldChange_table_2 %>% filter(!is.na(log2FoldChange) & !is.na(padj))



# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_2 <- log2FoldChange_table_2%>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05)

# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
diff_genes_2_immune <- diff_genes_2[diff_genes_2$geneName %in% hs_immune_genes$ensembl_gene, ]








########### 3 hours ############
design_matrix_3 <- filter(design_matrix, Timepoint %in% c(3))

count_matrix_3 <- count_matrix1[colnames(count_matrix1) %in% rownames(design_matrix_3)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_3) %in% rownames(design_matrix_3))

# Reorder the columns of count_matrix_3 to match the row names of design_matrix_3
count_matrix_3 <- count_matrix_3[, rownames(design_matrix_3)]
# order check
all(colnames(count_matrix_3) == rownames(design_matrix_3))

# Converting all read counts into integers
count_matrix_3[] <- lapply(count_matrix_3, as.integer)

design_matrix_3$variant <- factor(design_matrix_3$variant, levels = c("mock","delta","omicron"))
design_matrix_3$batch <- factor(design_matrix_3$batch)

dds_variant_time_3 <- DESeqDataSetFromMatrix(countData = count_matrix_3,
                                             colData = design_matrix_3,
                                             design = ~ batch + variant)

# set the factor level
dds_variant_time_3$variant <- relevel(dds_variant_time_3$variant, ref = "mock")
dds_variant_time_3$variant


# DESeq on 1h
deseq_variant_time_3 <- DESeq(dds_variant_time_3)


# Get differential expression results
result_delta_3 <- results(deseq_variant_time_3, contrast = c("variant","delta","mock"))
summary(result_delta_3)
head(result_delta_3)
result_delta_3 <- as.data.frame(result_delta_3)
result_delta_3 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_omicron_3 <- results(deseq_variant_time_3, contrast = c("variant","omicron","mock"))
summary(result_omicron_3)
head(result_omicron_3)
result_omicron_3 <- as.data.frame(result_omicron_3)
result_omicron_3 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_delta_3 <- result_delta_3[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_3")
df_omicron_3 <- result_omicron_3[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_3")

log2FoldChange_table_3 <- rbind(df_delta_3, df_omicron_3)

log2FoldChange_table_3 <- log2FoldChange_table_3 %>% filter(!is.na(log2FoldChange) & !is.na(padj))



# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_3 <- log2FoldChange_table_3%>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05)

# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
diff_genes_3_immune <- diff_genes_3[diff_genes_3$geneName %in% hs_immune_genes$ensembl_gene, ]






######################################################################################################################
############ Volcano ###########################################################################################

### 1 dpi
log2FoldChange_table_1$group <- factor(log2FoldChange_table_1$group, levels = c("Delta_1","Omicron_1"))

sig_alpha <- 0.05

day1 <- ggplot2::ggplot(log2FoldChange_table_1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5, na.rm = TRUE,
             aes(color = ifelse(log2FoldChange > 2 & -log10(padj) > -log10(sig_alpha), "red",
                                ifelse(log2FoldChange < -2 & -log10(padj) > -log10(sig_alpha), "dodgerblue", "darkgrey")))) +
  scale_color_manual(values = c("red" = "red", "dodgerblue" = "dodgerblue", "darkgrey" = "darkgrey"),
                     breaks = c("red", "dodgerblue"),
                     labels = c("Upregulated", "Downregulated"),
                     name = "Padj < 0.05 & FC > |2|") +
  theme_minimal() +
  xlim(c(-15, 15)) +
  ylim(c(-1, 10)) +
  ggtitle("log2 Fold Change for host genes in ALI cultures infected with SARS-CoV-2 VOCs compared to healthy controls at:
          1 dpi") +
  xlab("Log2 Fold Change") +
  ylab("-Log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5, size =14, face ="bold"),
        axis.text = element_text(size = 12, face = "italic"), 
        axis.title = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),  # Add panel border
        strip.background = element_blank(),  # Remove strip background
        panel.spacing = unit(3, "lines"),  # Reduce panel spacing
        legend.position = "bottom") +  
  facet_wrap(~ group, nrow = 1)


### 2 dpi
log2FoldChange_table_2$group <- factor(log2FoldChange_table_2$group, levels = c("Delta_2","Omicron_2"))

day2 <- ggplot2::ggplot(log2FoldChange_table_2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5, na.rm = TRUE,
             aes(color = ifelse(log2FoldChange > 2 & -log10(padj) > -log10(sig_alpha), "red",
                                ifelse(log2FoldChange < -2 & -log10(padj) > -log10(sig_alpha), "dodgerblue", "darkgrey")))) +
  scale_color_manual(values = c("red" = "red", "dodgerblue" = "dodgerblue", "darkgrey" = "darkgrey"),
                     breaks = c("red", "dodgerblue"),
                     labels = c("Upregulated", "Downregulated"),
                     name = "Padj < 0.05 & FC > |2|") +
  theme_minimal() +
  xlim(c(-15, 15)) +
  ylim(c(-1, 10)) +
  ggtitle("2 dpi") +
  xlab("Log2 Fold Change") +
  ylab("-Log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5, size =14, face ="bold"),
        axis.text = element_text(size = 12, face = "italic"), 
        axis.title = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),  # Add panel border
        strip.background = element_blank(),  # Remove strip background
        panel.spacing = unit(3, "lines"),  # Reduce panel spacing
        legend.position = "bottom") +  
  facet_wrap(~ group, nrow = 1)



#### 3 dpi
log2FoldChange_table_3$group <- factor(log2FoldChange_table_3$group, levels = c("Delta_3","Omicron_3"))

day3 <- ggplot2::ggplot(log2FoldChange_table_3, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5, na.rm = TRUE,
             aes(color = ifelse(log2FoldChange > 2 & -log10(padj) > -log10(sig_alpha), "red",
                                ifelse(log2FoldChange < -2 & -log10(padj) > -log10(sig_alpha), "dodgerblue", "darkgrey")))) +
  scale_color_manual(values = c("red" = "red", "dodgerblue" = "dodgerblue", "darkgrey" = "darkgrey"),
                     breaks = c("red", "dodgerblue"),
                     labels = c("Upregulated", "Downregulated"),
                     name = "Padj < 0.05 & FC > |2|") +
  theme_minimal() +
  xlim(c(-15, 15)) +
  ylim(c(-1, 10)) +
  ggtitle("3 dpi") +
  xlab("Log2 Fold Change") +
  ylab("-Log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5, size =14, face ="bold"),
        axis.text = element_text(size = 12, face = "italic"), 
        axis.title = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),  # Add panel border
        strip.background = element_blank(),  # Remove strip background
        panel.spacing = unit(3, "lines"),  # Reduce panel spacing
        legend.position = "bottom") +  
  facet_wrap(~ group, nrow = 1)



volcano_dpi <-  grid.arrange(day1 + theme(legend.position = "none"),
                             day2 + theme(legend.position = "none"),
                             day3 + theme(legend.position = "bottom"),
                             nrow = 3, heights = c(1,1,1.1))



######################################################################################################################



######################################################################################################################
######################################################################################################################
    ####################### GSEA for top up and downregulated genes for each VOC  ###############################
######################################################################################################################
######################################################################################################################
diff_genes_1 <- diff_genes_1 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_1$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))

diff_genes_1_immune <- diff_genes_1_immune %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_1_immune$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))
#######################################
# Beginning of infection: 
#######################################
### Omicron ###
diff_genes_1_Omicron <- diff_genes_1 %>%
  filter(group == "Omicron_1") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_1_Delta <- diff_genes_1 %>%
  filter(group == "Delta_1") %>%
  filter(!is.na(ENTREZID)) 



# Immune only
### Omicron ###
diff_genes_1_Omicron_immune <- diff_genes_1_immune %>%
  filter(group == "Omicron_1") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_1_Delta_immune <- diff_genes_1_immune %>%
  filter(group == "Delta_1") %>%
  filter(!is.na(ENTREZID)) 




#############
### Delta ###
# Ranking genes for GSEA
result_delta_1 <- result_delta_1 %>% filter(!is.na(log2FoldChange))
result_delta_1 <- result_delta_1 %>%
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_delta_1$geneName,
                           fromType = "ENSEMBL", toType = "ENTREZID",
                           column = "ENTREZID", keytype = "ENSEMBL"))

gene_list_d_1<- diff_genes_1_Delta_immune$log2FoldChange
names(gene_list_d_1) <- diff_genes_1_Delta_immune$ENTREZID
gene_list_d_1 = sort(gene_list_d_1, decreasing = TRUE)


#############
### Omicron ###
# Ranking genes for GSEA
result_omicron_1 <- result_omicron_1 %>% filter(!is.na(log2FoldChange))
result_omicron_1 <- result_omicron_1 %>%
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_omicron_1$geneName,
                           fromType = "ENSEMBL", toType = "ENTREZID",
                           column = "ENTREZID", keytype = "ENSEMBL"))

gene_list_o_1<- diff_genes_1_Omicron_immune$log2FoldChange
names(gene_list_o_1) <- diff_genes_1_Omicron_immune$ENTREZID
gene_list_o_1= sort(gene_list_o_1, decreasing = TRUE)







# ClusterprofilR GSEA GO pathway
gse_d_1 <- clusterProfiler::GSEA(gene_list_d_1,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               minGSSize = 15,
                               gson = NULL,
                               TERM2GENE = hs_gene2name_GO)


GOplot_d_1 <- ggplot(data = gse_d_1@result%>%
                     top_n(15, rank), 
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
  ggtitle("Delta 1 DPI") +
  scale_y_discrete(labels = custom_labeller)
GOplot_d_1




# ClusterprofilR GSEA REACTOME pathway
gse_o_1 <- clusterProfiler::GSEA(gene_list_o_1,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               minGSSize = 15,
                               gson = NULL,
                               TERM2GENE = hs_gene2name_GO)


GOplot_o_1 <- ggplot(data = gse_o_1@result%>%
                     top_n(15, rank), 
                   mapping = aes(x = NES, y = Description, 
                                 color = pvalue, size = setSize)) +
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
  ggtitle("Omicron 1 DPI") +
  scale_y_discrete(labels = custom_labeller)

GOplot_o_1



GO_plot_1 <-  grid.arrange(GOplot_d_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
                           GOplot_o_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "right"),
                           nrow = 1, ncol = 2, widths = c(1.1,1.15))


gsea_results_d_1 <- as.data.frame(gse_d_1@result)
gsea_results_d_1 <- as.data.frame(cbind(ID=gsea_results_d_1$ID,
                                        NES=gsea_results_d_1$NES,
                                        p.adjust=gsea_results_d_1$p.adjust))

gsea_results_o_1 <- as.data.frame(gse_o_1@result)
gsea_results_o_1 <- as.data.frame(cbind(ID=gsea_results_o_1$ID,
                                        NES=gsea_results_o_1$NES,
                                        p.adjust=gsea_results_o_1$p.adjust))


######################################
# 2nd day of infection:
######################################

diff_genes_2 <- diff_genes_2 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_2$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))

diff_genes_2_immune <- diff_genes_2_immune %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_2_immune$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))
#######################################
# Beginning of infection: 
#######################################
### Omicron ###
diff_genes_2_Omicron <- diff_genes_2 %>%
  filter(group == "Omicron_2") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_2_Delta <- diff_genes_2 %>%
  filter(group == "Delta_2") %>%
  filter(!is.na(ENTREZID)) 


# Immune only
### Omicron ###
diff_genes_2_Omicron_immune <- diff_genes_2_immune %>%
  filter(group == "Omicron_2") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_2_Delta_immune <- diff_genes_2_immune %>%
  filter(group == "Delta_2") %>%
  filter(!is.na(ENTREZID)) 


#############
### Delta ###
# Ranking genes for GSEA
result_delta_2 <- result_delta_2 %>% filter(!is.na(log2FoldChange))
result_delta_2 <- result_delta_2 %>%
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_delta_2$geneName,
                           fromType = "ENSEMBL", toType = "ENTREZID",
                           column = "ENTREZID", keytype = "ENSEMBL"))
gene_list_d_2<- result_delta_2$log2FoldChange
names(gene_list_d_2) <- result_delta_2$ENTREZID
gene_list_d_2= sort(gene_list_d_2, decreasing = TRUE)


#############
### Omicron ###
# Ranking genes for GSEA
result_omicron_2 <- result_omicron_2 %>% filter(!is.na(log2FoldChange))
result_omicron_2 <- result_omicron_2 %>%
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_omicron_2$geneName,
                           fromType = "ENSEMBL", toType = "ENTREZID",
                           column = "ENTREZID", keytype = "ENSEMBL"))
gene_list_o_2<- result_omicron_2$log2FoldChange
names(gene_list_o_2) <- result_omicron_2$ENTREZID
gene_list_o_2= sort(gene_list_o_2, decreasing = TRUE)




# ClusterprofilR GSEA REACTOME pathway
gse_d_2 <- clusterProfiler::GSEA(gene_list_d_2,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 minGSSize = 15,
                                 gson = NULL,
                                 TERM2GENE = hs_gene2name_GO)


GOplot_d_2 <- ggplot(data = gse_d_2@result%>%
                       top_n(15, NES), 
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
  ggtitle("Delta 2 DPI") +
  scale_y_discrete(labels = custom_labeller)
GOplot_d_2





# ClusterprofilR GSEA REACTOME pathway
gse_o_2 <- clusterProfiler::GSEA(gene_list_o_2,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 minGSSize = 15,
                                 gson = NULL,
                                 TERM2GENE = hs_gene2name_GO)


GOplot_o_2 <- ggplot(data = gse_o_2@result%>%
                       top_n(15, NES), 
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
  ggtitle("Omicron 2 DPI") +
  scale_y_discrete(labels = custom_labeller)

GOplot_o_2



GO_plot_2 <-  grid.arrange(GOplot_d_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
                           GOplot_o_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "right"),
                           nrow = 1, ncol = 2, widths = c(1.1,1.15))



gsea_results_d_2 <- as.data.frame(gse_d_2@result)
gsea_results_d_2 <- as.data.frame(cbind(ID=gsea_results_d_2$ID,
                                        NES=gsea_results_d_2$NES,
                                        p.adjust=gsea_results_d_2$p.adjust))

gsea_results_o_2 <- as.data.frame(gse_o_2@result)
gsea_results_o_2 <- as.data.frame(cbind(ID=gsea_results_o_2$ID,
                                        NES=gsea_results_o_2$NES,
                                        p.adjust=gsea_results_o_2$p.adjust))



######################################
# 3rd day of infection:
######################################

diff_genes_3 <- diff_genes_3 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_3$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))

diff_genes_3_immune <- diff_genes_3_immune %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_3_immune$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))
#######################################
# Beginning of infection: 
#######################################
### Omicron ###
diff_genes_3_Omicron <- diff_genes_3 %>%
  filter(group == "Omicron_3") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_3_Delta <- diff_genes_3 %>%
  filter(group == "Delta_3") %>%
  filter(!is.na(ENTREZID)) 

# Immune only
### Omicron ###
diff_genes_3_Omicron_immune <- diff_genes_3_immune %>%
  filter(group == "Omicron_3") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_3_Delta_immune <- diff_genes_3_immune %>%
  filter(group == "Delta_3") %>%
  filter(!is.na(ENTREZID)) 




#############
### Delta ###
# Ranking genes for GSEA
result_delta_3 <- result_delta_3 %>% filter(!is.na(log2FoldChange))
result_delta_3 <- result_delta_3 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_delta_3$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))

gene_list_d_3<- diff_genes_3_Delta_immune$log2FoldChange
names(gene_list_d_3) <- diff_genes_3_Delta_immune$ENTREZID
gene_list_d_3= sort(gene_list_d_3, decreasing = TRUE)


#############
### Omicron ###
# Ranking genes for GSEA
result_omicron_3 <- result_omicron_3 %>% filter(!is.na(log2FoldChange))
result_omicron_3 <- result_omicron_3 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_omicron_3$geneName, 
                           fromType = "ENSEMBL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "ENSEMBL"))

gene_list_o_3<- diff_genes_3_Omicron_immune$log2FoldChange
names(gene_list_o_3) <- diff_genes_3_Omicron_immune$ENTREZID
gene_list_o_3= sort(gene_list_o_3, decreasing = TRUE)





# ClusterprofilR GSEA REACTOME pathway
gse_d_3 <- clusterProfiler::GSEA(gene_list_d_3,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 minGSSize = 15,
                                 gson = NULL,
                                 TERM2GENE = hs_gene2name_GO)


GOplot_d_3 <- ggplot(data = gse_d_3@result%>%
                       top_n(15, setSize), 
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
  ggtitle("Delta 3 DPI") +
  scale_y_discrete(labels = custom_labeller)
GOplot_d_3






# ClusterprofilR GSEA REACTOME pathway
gse_o_3 <- clusterProfiler::GSEA(gene_list_o_3,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 minGSSize = 15,
                                 gson = NULL,
                                 TERM2GENE = hs_gene2name_GO)


GOplot_o_3 <- ggplot(data = gse_o_3@result%>%
                       top_n(15, NES), 
                     mapping = aes(x = rank, y = Description, 
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
  ggtitle("Omicron 3 DPI") +
  scale_y_discrete(labels = custom_labeller)

GOplot_o_3



GO_plot_3 <-  grid.arrange(GOplot_d_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
                           GOplot_o_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "right"),
                           nrow = 1, ncol = 2, widths = c(1.1,1.15))



gsea_results_d_3 <- as.data.frame(gse_d_3@result)
gsea_results_d_3 <- as.data.frame(cbind(ID=gsea_results_d_3$ID,
                                        NES=gsea_results_d_3$NES,
                                        p.adjust=gsea_results_d_3$p.adjust))

gsea_results_o_3 <- as.data.frame(gse_o_3@result)
gsea_results_o_3 <- as.data.frame(cbind(ID=gsea_results_o_3$ID,
                                        NES=gsea_results_o_3$NES,
                                        p.adjust=gsea_results_o_3$p.adjust))

grid.arrange(
  GO_plot_1,
  GO_plot_2,
  GO_plot_3 ,  
  ncol = 1, nrow= 3)






