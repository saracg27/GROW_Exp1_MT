library(EBSeq)
library(here)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)


# Set the directory where your files are located
data_dir <- here("data", "contrasts") 
# Get a list of all the files in the directory
file_list <- list.files(path = data_dir, pattern = "^merged_gene_abundance_D0_contrast.*\\.rds$", full.names = TRUE) 
# Get the total number of files
total_files <- length(file_list)


# Iterate over each file
for (file in file_list) {
   
   # Extract the file name without the extension
   file_name <- tools::file_path_sans_ext(basename(file))
   # Extract the number from the file name
   file_number <- str_extract(file_name, "\\d+")
   # Construct the metadata file name
   meta_file <- file.path(data_dir, paste0("metaD0_contrast", file_number, ".rds"))
   # Extract the file name without the extension
   meta_name <- tools::file_path_sans_ext(basename(meta_file))
   # Check if the corresponding metadata file exists
   if (!file.exists(meta_file)) {
      print("Metadata file not found for", file, "\n")
      next
   }
   
   
   # Load the data from the rds object
   #countData <- vroom(file = file, col_names = TRUE)
   countData <- readRDS(file = file)
   countData <- countData[,order(colnames(countData))]

   # Load the metadata from the csv file
   meta <- readRDS(file = meta_file)
   meta <- meta[order(rownames(meta)),]

      # Sanity Check 
   if (all(row.names(meta) == colnames(countData))) {
      print("Sanity check ok: meta rownames same as countData colnames")
   } else {
      print("ERROR: not all row names are the same")
   }
   
   #
   
   
   # Create and subset splitting of contrasts
   for (i in seq_along(contrast)) {
      # Subset by contrast
      meta_i <- meta[meta$treatment %in% treatments[[i]], ]
      common_names <- intersect(colnames(countData), rownames(meta_i))
      countData_i <- countData[, common_names]
      # Filter out zeros
      countData_i <- countData_i[!rowSums(countData_i) == 0, ]
      # Count the number of rows
      row_count_i <- nrow(countData_i)
      dim(row_count_i)
      # Print the row count
      #cat("The contrast ", countData_i, "from", file_name, " has ", row_count_i, "rows")
      
      # Save the data as RDS, tsv, and csv files
      saveRDS(countData_i, file = here("data", "clean", paste0(file_name, "_", contrast[i], "_0out.rds")))
      saveRDS(meta_i, file = here("data", "clean", paste0(meta_name, "_", contrast[i], "_0out.rds")))
      write.table(countData_i, file = here("data", "clean", paste0(file_name, "_", contrast[i], "_0out.tsv")), sep = "\t", quote = FALSE)
      write.csv(meta_i, file = here("data", "clean", paste0(meta_name, "_", contrast[i], "_0out.csv")), row.names = FALSE)
      
      # Remove the heavy tsv file and its corresponding metadata from the environment
      rm(countData_i, meta_i)
      
      # Print a message indicating the completion of processing for the current contrast
      #cat("Subsetting completed for", contrast[i], "out of", length(contrast), "-", file, "\n")
   }
   
   
   #Save memory
   rm(countData)
   rm(meta)
   
   # Print a message indicating the completion of processing for the current file
   #print("Processing completed for", file, "\n")
}











#Data file should be genes in rows and samples in columns

#EB test from EB seq for multiple testing of hypothesis

#rhizo 25 vs.100
map.s.rhizo <- map.s[(map.s$Sample == "Rhizosphere") & (map.s$WHC == "25" | map.s$WHC == "100"), ] # 12 obs in 3 var
MT.rhizo <- MT.s[(map.s$Sample == "Rhizosphere") & (map.s$WHC == "25" | map.s$WHC == "100"), ]


MT.rhizo <- t(MT.rhizo)


dim(MT.rhizo) #

#Create group file L = 25, H = 100
SWC.rhizo = as.factor(c("H", "L","H", "L","H", "L","H", "L", "L", "H", "L", "H"))

#Create size Factors file
size.rhizo = MedianNorm(MT.rhizo)

#Run EBTest
EBtestOUT.rhizo <- EBTest(MT.rhizo, Conditions = SWC.rhizo, sizeFactors = size.rhizo, maxround = 5)

#Get results
DE.rhizo <- GetDEResults(EBtestOUT.rhizo, FDR = 0.05)

#Calculate log2 fold change and mean counts
fold.rhizo <- data.frame(row.names = row.names(MT.rhizo))
fold.rhizo$mean <- rowMeans(MT.rhizo)
fold.rhizo$log2 <- log2(rowMeans(MT.rhizo[,SWC.rhizo == "L"])/rowMeans(MT.rhizo[,SWC.rhizo == "H"]))

#Get annotations for DE at 0.05
rhizo.DE.annot <- annot[annot$gene_id %in% DE.rhizo$DEfound,]

#Get annotations for DE
DE.bact.rhizo <- rhizo.DE.annot[rhizo.DE.annot$tax_kingdom == "Bacteria",2]
bact.point.rhizo <- fold.rhizo[row.names(fold.rhizo) %in% DE.bact.rhizo,]
DE.arch.rhizo <- rhizo.DE.annot[rhizo.DE.annot$tax_kingdom == "Archaea",2]
arch.point.rhizo <- fold.rhizo[row.names(fold.rhizo) %in% DE.arch.rhizo,]
DE.fun.rhizo <- rhizo.DE.annot[grep("mycota", rhizo.DE.annot$tax_phylum, ignore.case = TRUE),2]
fun.point.rhizo <- fold.rhizo[row.names(fold.rhizo) %in% DE.fun.rhizo,]
DE.plant.rhizo <- rhizo.DE.annot[rhizo.DE.annot$tax_phylum == "Streptophyta",2]
plant.point.rhizo <- fold.rhizo[row.names(fold.rhizo) %in% DE.plant.rhizo,]
DE.null.rhizo <- rhizo.DE.annot[rhizo.DE.annot$tax_kingdom == "NULL",2]
null.point.rhizo <- fold.rhizo[row.names(fold.rhizo) %in% DE.null.rhizo,]
DE.other.rhizo <- rhizo.DE.annot[!rhizo.DE.annot$gene_id %in% append(append(DE.bact.rhizo,DE.fun.rhizo),append(DE.plant.rhizo, append(DE.null.rhizo,DE.arch.rhizo))),2]
other.point.rhizo <- fold.rhizo[row.names(fold.rhizo) %in% DE.other.rhizo,]

#Plots
#Only plot first 100 000 points for clarity
volcano.rhizo <- ggplot() +
   geom_point(data = fold.rhizo[1:100000,], aes(x = mean, y = log2)) +
   geom_point(data = other.point.rhizo, aes(x = mean, y = log2), color = "yellow") +
   geom_point(data = null.point.rhizo, aes(x = mean, y = log2), color = "purple") +
   geom_point(data = bact.point.rhizo, aes(x = mean, y = log2), color = "red") +
   geom_point(data = fun.point.rhizo, aes(x = mean, y = log2), color = "blue") +
   geom_point(data = arch.point.rhizo, aes(x = mean, y = log2), color = "pink") +
   geom_point(data = plant.point.rhizo, aes(x = mean, y = log2), color = "green") +
   geom_hline(yintercept = 0, linetype = "solid") +
   xlab("Mean read count") + 
   ylab("Log 2 fold change (25%/100%)") + 
   scale_x_log10() +
   theme_bw()
volcano.rhizo 

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
