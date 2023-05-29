
### to work on a cluster, need to load packages within script, installed using R interactive mode
# The code combines data loading, manipulation, and file I/O operations in a loop over multiple files. 
# It subsets the data based on contrasts and filters out rows with zero counts. It also performs some data cleaning 
# and saves the resulting data and metadata in various formats. Previous step to DE
### Run on the cluster using split_contrasts.sh from the code folder: 

## 1) Navigate on globus to the desired folder : /home/saracoga/projects/def-yergeaue-ab/saracoga/GROW project/mesocosm1_MT/
## 2) transfer both the Split_contrasts.R and the split_contrasts.sh scripts
#### In the terminal: 
## 3) connect to graham on your terminal: ssh -vvv -Y saracoga@graham.computecanada.ca
## 4) Navigate to the folder: cd projects/def-yergeaue-ab/saracoga/GROW\ project/mesocosm1_MT/
## 5) Submit the job: sbatch split_contrasts.sh

### Info for first submission: 
#sbatch: NOTE: Your memory request of 131072M was likely submitted as 128G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
#Submitted batch job 6811168




## sbatch split_contrasts.sh


library(data.table)
library(here)
library(vroom)
library(tidyverse)

# Set the contrast and their corresponding treatments
contrast <- c("contrast1", "contrast2", "contrast3", "contrast4", "contrast5", "contrast6")
treatments <- list(c("sdmt_unplanted_ctl", "sdmt_unplanted_ospw"),
                   c("sdmt_unplanted_ctl", "sdmt_carex_ospw"),
                   c("sdmt_unplanted_ospw", "sdmt_carex_ospw"),
                   c("sdmt_carex_ospw", "rhizo_carex_ospw"),
                   c("sdmt_unplanted_ospw", "rhizo_carex_ospw"),
                   c("sdmt_unplanted_ctl", "rhizo_carex_ospw"))
treatments1 <- c("sdmt_unplanted_ctl - sdmt_unplanted_ospw", "sdmt_unplanted_ctl - sdmt_carex_ospw", "sdmt_unplanted_ospw - sdmt_carex_ospw", "sdmt_carex_ospw - rhizo_carex_ospw", "sdmt_unplanted_ospw - rhizo_carex_ospw", "sdmt_unplanted_ctl - rhizo_carex_ospw")
contrast_table <- as.data.frame(cbind(contrast, treatments1))
write.csv(contrast_table, file = here("data", "clean", "contrast_checked.csv"), row.names = FALSE)
rm(contrast_table)
rm(treatments1)


# Set the directory where your files are located
data_dir <- here("data", "clean") 
# Get a list of all the files in the directory
file_list <- list.files(path = data_dir, pattern = "^merged_gene_abundance_D.*\\.tsv$", full.names = TRUE) 
# Get the total number of files
total_files <- length(file_list)


# Iterate over each file
for (file in file_list) {
    
    # Extract the file name without the extension
    file_name <- tools::file_path_sans_ext(basename(file))
    # Extract the number from the file name
    file_number <- str_extract(file_name, "\\d+")
    # Construct the metadata file name
    meta_file <- file.path(data_dir, paste0("metaD", file_number, ".csv"))
    # Extract the file name without the extension
    meta_name <- tools::file_path_sans_ext(basename(meta_file))
    # Check if the corresponding metadata file exists
    if (!file.exists(meta_file)) {
        print("Metadata file not found for", file, "\n")
        next
    }

    
    # Load the data from the tsv file
    countData <- vroom(file = file, col_names = TRUE)
    # Clean countData names to match metadata
    countData <- as.data.frame(countData) # As dataframe to reoder in base R
    colnames(countData) <- sub(".*_i5_", "", colnames(countData))
    colnames(countData) <- sub("^[^.]+\\.", "", colnames(countData))
    countData <- countData[,order(colnames(countData))]
    # Change row names to feature_id
    row.names(countData) <- countData$feature_id
    countData$feature_id <- NULL
    
    
    # Load the metadata from the csv file
    meta <- read.csv(meta_file)
    # Keep only samples present in countData
    meta <- meta[!meta$Compartment == "Roots",]
    meta$Time <- NULL
    meta$Sampling_date <- NULL
    meta$SampleID <- NULL
    meta$Compartment <- factor(meta$Compartment, levels = c("Sediments", "Rhizosphere"), labels = c("sdmt", "rhizo"))
    meta$Sample_type <- factor(meta$Sample_type,levels = c("No_plant", "Carex"), labels = c("unplanted", "carex"))
    meta$Water_type <- factor(meta$Water_type, levels = c("Artificial_OSPW", "OSPW"), labels = c("ctl", "ospw"))
    meta$treatment <- paste0(meta$Compartment, "_", meta$Sample_type, "_", meta$Water_type)
    meta$treatment <- factor(meta$treatment,  levels = c("sdmt_unplanted_ctl", "sdmt_unplanted_ospw", "sdmt_carex_ospw", "rhizo_carex_ospw" ))
    meta <- meta[order(meta$alias),]
    rownames(meta) <- meta$alias
    
    # Sanity Check 
    if (all(row.names(meta) == colnames(countData))) {
        print("Sanity check ok: meta rownames same as countData colnames")
    } else {
        print("ERROR: not all row names are the same")
    }
    

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

