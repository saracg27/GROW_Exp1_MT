### DE analysis for D-0
### Run on the cluster using de_d0.sh from the code folder: sbatch de_d0.sh
### to work on a cluster, need to load packages within script, installed using R interactive mode

library(vroom)
library(data.table)
library(here)
library(tidyverse)
library(gdata)
library(DESeq2)
# Data
countData <- vroom(file = here::here('data', 'clean', 'merged_gene_abundance_D-0.tsv'), col_names = TRUE)
meta <- read.csv(file = here::here('data', 'clean', 'metaD0.csv'), header = TRUE, sep = ",")

# Clean Data

## Clean countData names to match metadata
countData <- as.data.frame(countData) # As dataframe to reoder in base R
colnames(countData) <- sub(".*_i5_", "", colnames(countData))
colnames(countData) <- sub("^[^.]+\\.", "", colnames(countData))
countData <- countData[,order(colnames(countData))]
names(countData)
## Change row names to feature_id
row.names(countData) <- countData$feature_id
countData$feature_id <- NULL
## Keep only samples present in countData (not root samples)
meta <- meta[!meta$Compartment == "Roots",]
## Recode factors in meta
meta$Time <- NULL
meta$Sampling_date <- NULL
meta$SampleID <- NULL
meta$Compartment <- factor(meta$Compartment,
                           levels = c("Sediments", "Rhizosphere"),
                           labels = c("sdmt", "rhizo"))
meta$Sample_type <- factor(meta$Sample_type,
                           levels = c("No_plant", "Carex"),
                           labels = c("unplanted", "carex"))
meta$Water_type <- factor(meta$Water_type,
                          levels = c("Artificial_OSPW", "OSPW"),
                          labels = c("ctl", "ospw"))
meta$treatment <- paste0(meta$Compartment, "_", meta$Sample_type, "_", meta$Water_type)
meta$treatment <- factor(meta$treatment, 
                         levels = c("sdmt_unplanted_ctl", "sdmt_unplanted_ospw", "sdmt_carex_ospw", "rhizo_carex_ospw" ))
meta <- meta[order(meta$alias),]
rownames(meta) <- meta$alias

# Inspect Data

## Summary
summary <- summary(countData) 
par(mar = c(15,5,4,1) + 2)
library_size <- colSums(countData)
library_size_plot <- barplot(colSums(countData)/1e6, las = 3)
#histogram_raw <- hist(countData$`25-Mesocosm-1-1-Exp1-D8-Sed`, br = 100)

## Log2 transform and visualize
logCountData = log2(1 + countData)
#histogram_log <- hist(logCountData$`25-Mesocosm-1-1-Exp1-D8-Sed`, br = 50)
plot_log <- plot(logCountData[,1], logCountData[,2])

## Reduce uninformative features Threashold of 10 (decided somewhat randomly)
countData1 <- countData[rowSums(countData) > 10,] # Down from 26M genes to 2 163 717  genes (2.16M)

## Check data is in order

if (all(row.names(meta) == colnames(countData1))) {
   print("Sanity check ok")
} else {
   print("ERROR: not all row names are equal to col names same")
}

# Create Deseq object
## Create the DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = countData1,
                              colData = meta, 
                              design = ~ treatment)

## Create the DESeq object
dds <- DESeq(dds) # Takes 90 minutes approximately
nrow(dds)

## Look at size factors
size_factors <- sizeFactors(dds)

## Apply a regularized log transformation
rld <- rlog(dds)
plotpca <- plotPCA(rld, intgroup = "treatment")

save.image(file = here("output", "DE_D0.RData"))
saveRDS(dds, file = here("output", "dds_d0.rds" ))
           