library(EBSeq)
library(here)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

meta <- readRDS(file = here('data', 'contrasts', 'metaD0_contrast1.rds'))
#Data file should be genes in rows and samples in columns - Rows with no zeros
countData <- readRDS(file = here('data', 'contrasts', 'merged_gene_abundance_D0_contrast1.rds'))


# Sanity check for Zero filter and order of columns and rownames
if (all(rowSums(countData) != 0)) {
   print("rowSums zero filter check ok")
} else {
   print("ERROR: there are rowSums == 0 - Need to filter out non significant features")
}
length(which(rowSums(countData) == 0)) # Should be zero
# Sanity Check 
if (all(row.names(meta) == colnames(countData))) {
   print("Sanity check ok: meta rownames same as countData colnames")
} else {
   print("ERROR: not all row names are the same")
}
countData2 <- countData1[rowSums(countData1) > 25000,] 

# Create group file
treatment <- as.vector(meta$treatment)
treatment
treat_factor <- as.factor(paste0(treatment, sep = "|"))
treat_factor <- as.factor(c("sdmt_unplanted_ospw", "sdmt_unplanted_ctl", "sdmt_unplanted_ctl",  "sdmt_unplanted_ctl",  "sdmt_unplanted_ospw", "sdmt_unplanted_ospw",
                            "sdmt_unplanted_ospw" ,"sdmt_unplanted_ospw", "sdmt_unplanted_ospw" ,"sdmt_unplanted_ctl" , "sdmt_unplanted_ctl" , "sdmt_unplanted_ospw",
                            "sdmt_unplanted_ospw", "sdmt_unplanted_ctl" , "sdmt_unplanted_ctl" , "sdmt_unplanted_ctl" ))
# Create size Factors
size_treatment <- MedianNorm(countData2)
countData2 <- as.matrix(countData2)
#Run EBTest
EBtestOUT <- EBTest(countData2, Conditions = treat_factor, sizeFactors = size_treatment, maxround = 7)

#Get results
DE <- GetDEResults(EBtestOUT, FDR = 0.05)

#Calculate log2 fold change and mean counts
fold <- data.frame(row.names = row.names(countData2))
fold$mean <- rowMeans(countData2)
fold$log2 <- log2(rowMeans(countData2[,treat_factor == "sdmt_unplanted_ospw"])/rowMeans(countData2[,treat_factor == "sdmt_unplanted_ctl"]))

#Get annotations for DE at 0.05
annot = vroom(file = here::here("data", "raw", "annotations.tsv"))
annot2 = annot[annot$gene_id %in% DE$DEfound,]

#Get annotations for DE
DE.bact <- annot2[annot2$tax_kingdom == "k__Bacteria",2]
bact.point <- fold[row.names(fold) %in% DE.bact$gene_id,]
DE.metaz <- annot2[annot2$tax_kingdom == "k__Metazoa",2]
metaz.point <- fold[row.names(fold) %in% DE.metaz$gene_id,]
DE.fun <- annot2[grep("mycota", annot2$tax_phylum, ignore.case = TRUE),2]
fun.point <- fold[row.names(fold) %in% DE.fun$gene_id,]
DE.plant <- annot2[annot2$tax_kingdom == "k__Viridiplantae",2]
plant.point <- fold[row.names(fold) %in% DE.plant$gene_id,]
DE.null <- annot2[annot2$tax_kingdom == "k__NULL",2]
null.point <- fold[row.names(fold) %in% DE.null$gene_id,]
DE.other <- annot2[!annot2$gene_id %in% append(append(DE.bact,DE.fun),append(DE.plant, append(DE.null,DE.metaz))),2]
other.point <- fold[row.names(fold) %in% DE.other$gene_id,]


#Plots
#Only plot first 100 000 points for clarity
volcano <- ggplot() +
   geom_point(data = fold[1:100000,], aes(x = mean, y = log2)) +
   geom_point(data = null.point, aes(x = mean, y = log2), color = "purple") +
   geom_point(data = bact.point, aes(x = mean, y = log2), color = "red") +
   geom_point(data = fun.point, aes(x = mean, y = log2), color = "blue") +
   geom_point(data = metaz.point, aes(x = mean, y = log2), color = "yellow") +
   geom_point(data = plant.point, aes(x = mean, y = log2), color = "green") +
   geom_hline(yintercept = 0, linetype = "solid") +
   xlab("Mean read count") + 
   ylab("Log 2 fold change D0 sedmnt OSPW vs CTRL") + 
   scale_x_log10() +
   theme_bw()
volcano 

##### CORRECT FROM HERE ###########

#Sort for exporting
rhizo.DE.annot <- rhizo.DE.annot[order(rhizo.DE.annot$gene_id),]
DEfound.rhizo <- data.frame(DE.rhizo$DEfound)
DEfound.rhizo <- DEfound.rhizo[order(DEfound.rhizo$DE.rhizo.DEfound),]
fold.rhizo.DE <- fold.rhizo[row.names(fold.rhizo) %in% DE.rhizo$DEfound,]
fold.rhizo.DE <- fold.rhizo.DE[order(row.names(fold.rhizo.DE)),]
pval.rhizo.DE <- data.frame(DE.rhizo$PPMat[row.names(DE.rhizo$PPMat) %in% DE.rhizo$DEfound,])
pval.rhizo.DE <- pval.rhizo.DE[order(row.names(pval.rhizo.DE)),]
sum(rhizo.DE.annot$gene_id == DEfound.rhizo) #Sanity check, should be 21765
sum(row.names(fold.rhizo.DE) == DEfound.rhizo) #Sanity check, should be 21765
sum(row.names(pval.rhizo.DE) == DEfound.rhizo)#Sanity check, should be 21765
DEfound.rhizo <- cbind(DEfound.rhizo,fold.rhizo.DE, pval.rhizo.DE, rhizo.DE.annot)

#Export
#Matrix of all genes with DE and EE p-values
write.table(as.matrix(DE.rhizo$PPMat), file = here("output", "PPmat_rhizo.txt"), eol = "\n", sep = "\t")
#List of the DE transcripts, p-Values, annotations, abundance
write.table(DEfound.rhizo, file = here("output", "DEfound_rhizo.txt"), eol = "\n", sep = "\t")
#Status of the transcript: DE, EE, or Filtered+Reason
write.table(as.matrix(DE.rhizo$Status), file = here("output", "Status_rhizo.txt"), eol = "\n", sep = "\t")








