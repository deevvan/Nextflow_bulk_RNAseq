# script to perform WGCNA

setwd("/mmfs1/projects/changhui.yan/DeewanB/VOC_SRR/deseq2_human")

#devtools::install_github("kevinblighe/CorLevelPlot")


load(file = "WGCNA_data.RData")


library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Fetch Data ------------------------------------------------
#data <- read.delim('GSE152418_WGCNA_data.txt', header = T)
data <- read.csv("GSE157103_WGCNA_data.csv",header=T)

# get metadata
#geo_id <- "GSE152418"
geo_id <- "GSE157103"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
#phenoData_sample <- phenoData_sample[,c(1,2,46:50)]
phenoData <- phenoData[,c(1,2,10,13,37)]
#phenoData_original <- phenoData
#phenoData <- phenoData_original

# Corrected filtering command to remove specific rows
phenoData <- phenoData[!(phenoData$characteristics_ch1 == "disease state: non-COVID-19" &
                           phenoData$characteristics_ch1.3 == "icu: yes"), ]


# prepare data
data <- data %>% 
        gather(key = 'samples', value = 'counts', -X.symbol) %>% 
        inner_join(., phenoData, by = c('samples' = 'description')) %>% 
        dplyr::select(1,3,5) %>% 
        spread(key = 'geo_accession', value = 'counts') %>% 
        column_to_rownames(var = 'X.symbol') 


# 2. WGCNA's QC step for outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK # allOK = TRUE means all the SRR runs passed and there are no outliers

table(gsg$goodGenes) # check outliers in rows (genes)
table(gsg$goodSamples) # check outliers in columns. (SRR runs)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,] # Selecting only genes that are not outliers

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples
samples.to.be.excluded <- c('GSM4753087', 'GSM4753095', 'GSM4753139', 'GSM4753146')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# Rechecking outlier samples - hierarchical clustering - method 1 after excluding outliers
htree.subset <- hclust(dist(t(data.subset)), method = "average")
plot(htree.subset)

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)
names(colData) <- c("title","geo_accession","disease_state","ICU_state","description" )
#names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

data.subset[] <- sapply(data.subset, as.numeric)
data.subset <- round(data.subset)  # Round to nearest integer

#non_integer_values <- data.subset[!is.integer(data.subset)]
#any(!is.integer(data.subset))

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not specifying model



## remove all genes with counts < 15 in more than 75% of samples (122*0.75=91.5)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 91,]
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
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

# Determining Minimum mean connectivity 
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# Determined soft threshold 
soft_power <- 20 # at 20, Rsquare value was just above 80% line and mean connectivity was lowest consecutively

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

# Assigining WGCNA correlation to cor variable to avoid namespace error
temp_cor <- cor
cor <- WGCNA::cor

# Blockwise module functions's memory estimate with respect to blocksize (4GB:8k, 16GB:20k, 32GB:30k)

bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 30000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor # reassign the correlation funtion back the to original cor function 


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

# get number of genes for each module
table(bwnet$colors)

# Plot dendrogram and module colors before and after merging genes
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors), # visualize genes before and after merging
                    c("unmerged", "merged"),
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
          mutate(disease_state_bin = ifelse(!grepl('non', disease_state), 1, 0)) %>% 
          dplyr::select(7)


# Add a new column 'Severity' based on specified conditions
colData <- colData %>%
  mutate(Severity = ifelse(grepl('yes', ICU_state) & (!grepl('non', disease_state)), "severe", "normal"),
         .keep = "all")

# binarize categorical variable ICU status
severity <- colData %>%
  mutate(ICU_state_bin = ifelse(grepl('yes', ICU_state), 1, 0)) %>%
  dplyr::select(7)

# severity <- colData %>% 
#              mutate(ICU_state_bin = ifelse(grepl('severe', Severity), 1, 0)) %>% 
#              dplyr::select(7)

traits <- cbind(traits, severity)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)


# visualize module-trait association as a heatmap
# adding binary train classes to ME dataframe
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:19],
             y = names(heatmap.data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

module.gene.mapping <- as.data.frame(bwnet$colors)

# COPY OF MAIN MODULE MAP
module.gene.mapping1 <- module.gene.mapping %>%
  tibble::rownames_to_column(var = "RowNames")

ICU_EG_genelist <- module.gene.mapping %>% 
                    filter(`bwnet$colors` == 'black') %>% 
                    rownames() %>%
                    as.data.frame(.)

ICU_EG_genelist2 <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames() %>%
  as.data.frame(.)

disease_EG_genelist <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames() %>%
  as.data.frame(.)
write.csv(disease_EG_genelist,"turquoiseEM.csv")

# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and 
# the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


# Calculate the gene significance and associated p-values

gene.signf.corr.ICU <- cor(norm.counts, traits$ICU_state_bin, use = 'p')
gene.signf.corr.ICU.pvals <- corPvalueStudent(gene.signf.corr.ICU, nSamples)

# The top 25 genes based on pvalue that are associated with severity/ICU status
gene.signf.corr.ICU.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  #rownames(.) %>%
  head(25)



gene.signf.corr.disease <- cor(norm.counts, traits$disease_state_bin, use = 'p')
gene.signf.corr.disease.pvals <- corPvalueStudent(gene.signf.corr.disease, nSamples)

# The top 25 genes based on pvalue that are associated with severity/ICU status
gene.signf.corr.disease.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

########## Clustering based on severity alone in COVID-19 patients:

# Corrected filtering command to remove specific rows
phenoData <- phenoData[!(phenoData$characteristics_ch1 == "disease state: non-COVID-19"), ]


# prepare data
data <- data %>% 
  gather(key = 'samples', value = 'counts', -X.symbol) %>% 
  inner_join(., phenoData, by = c('samples' = 'description')) %>% 
  dplyr::select(1,3,5) %>% 
  spread(key = 'geo_accession', value = 'counts') %>% 
  column_to_rownames(var = 'X.symbol') 


# 2. WGCNA's QC step for outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK # allOK = TRUE means all the SRR runs passed and there are no outliers

table(gsg$goodGenes) # check outliers in rows (genes)
table(gsg$goodSamples) # check outliers in columns. (SRR runs)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,] # Selecting only genes that are not outliers

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples
samples.to.be.excluded <- c('GSM4753054', 'GSM4753070', 'GSM4753106', 'GSM4753093', 'GSM4753103','GSM4753108',
                            'GSM4753050','GSM4753100','GSM4753077','GSM4753081','GSM4753109','GSM4753061','GSM4753059',
                            'GSM4753068','GSM4753087','GSM4753095','GSM4753116','GSM4753120')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# Rechecking outlier samples - hierarchical clustering - method 1 after excluding outliers
htree.subset <- hclust(dist(t(data.subset)), method = "average")
plot(htree.subset)

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)
names(colData) <- c("title","geo_accession","disease_state","ICU_state","description" )
#names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

data.subset[] <- sapply(data.subset, as.numeric)
data.subset <- round(data.subset)  # Round to nearest integer

#non_integer_values <- data.subset[!is.integer(data.subset)]
#any(!is.integer(data.subset))

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (122*0.75=91.5)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 61,]
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
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

# Determining Minimum mean connectivity 
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# Determined soft threshold 
soft_power <- 20 # at 20, Rsquare value was just above 80% line and mean connectivity was lowest consecutively

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

# Assigining WGCNA correlation to cor variable to avoid namespace error
temp_cor <- cor
cor <- WGCNA::cor

# Blockwise module functions's memory estimate with respect to blocksize (4GB:8k, 16GB:20k, 32GB:30k)

bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 30000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor # reassign the correlation funtion back the to original cor function 


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

# get number of genes for each module
table(bwnet$colors)

# Plot dendrogram and module colors before and after merging genes
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors), # visualize genes before and after merging
                    c("unmerged", "merged"),
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
  mutate(ICU_state_bin = ifelse(grepl('yes', ICU_state), 1, 0)) %>% 
  dplyr::select(6)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)


# visualize module-trait association as a heatmap
# adding binary train classes to ME dataframe
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[16],
             y = names(heatmap.data)[1:15],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

module.gene.mapping <- as.data.frame(bwnet$colors)

# COPY OF MAIN MODULE MAP
module.gene.mapping1 <- module.gene.mapping %>%
  tibble::rownames_to_column(var = "RowNames")

ICU_EG_genelist <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'red') %>% 
  rownames() %>%
  as.data.frame(.)


#write.csv(ICU_EG_genelist,"redEM.csv")

# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and 
# the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


# Calculate the gene significance and associated p-values

gene.signf.corr.ICU <- cor(norm.counts, traits$ICU_state_bin, use = 'p')
gene.signf.corr.ICU.pvals <- corPvalueStudent(gene.signf.corr.ICU, nSamples)

# The top 25 genes based on pvalue that are associated with severity/ICU status
gene.signf.corr.ICU.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  #rownames(.) %>%
  head(25)



gene.signf.corr.disease <- cor(norm.counts, traits$disease_state_bin, use = 'p')
gene.signf.corr.disease.pvals <- corPvalueStudent(gene.signf.corr.disease, nSamples)

# The top 25 genes based on pvalue that are associated with severity/ICU status
gene.signf.corr.disease.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)






### Trial GS- MM plot

GS.trait = as.numeric(cor(norm.counts, traits, use = "p"))

moduleLabelsAutomatic = bwnet$colors
# Convert labels to colors for plotting
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
# Choose a module assignment
moduleColorsICU = moduleColorsAutomatic

# Define numbers of genes and samples
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(norm.counts, moduleColorsICU)$eigengenes
MEsICU = orderMEs(MEs0)

modTraitCor = cor(MEsICU, traits, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", 
                   sep = "")
dim(textMatrix) = dim(modTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = modTraitCor, xLabels = names(traits), yLabels = names(MEsICU), 
               ySymbols = names(MEsICU), colorLabels = FALSE, colors = greenWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), 
               main = paste("Module-trait relationships"))


# extract colnames of datKME from 4th character onwards (eg. kMEblack = black)
colorOfColumn = substring(names(datKME), 4)
par(mfrow = c(2, 2))

selectModules = c("greenyellow", "midnightblue", "black", "cyan")
par(mfrow = c(2, length(selectModules)/2))


datKME = signedKME(norm.counts, MEsICU)

for (module in selectModules) {
  # Finds corresponding column index in datKME that matches module name (or color)
  column = match(module, colorOfColumn)
  # Creates a logical vector indicating rows in moduleColorsICU that match the current module.
  restModule = moduleColorsICU == module
  # Subsets datKME to include rows corresponding to the current module and the specified column and plots
  verboseScatterplot(datKME[restModule, column], GS.trait[restModule], xlab = paste("Module Membership ", module, "module"), 
                     ylab = "GS", main = paste("kME.", module,"vs. GS"), col = module)
}


# Save the current R session
#save.image(file = "WGCNA_data.RData")




