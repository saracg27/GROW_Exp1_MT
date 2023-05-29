library(EBSeq)
library(here)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

meta <- readRDS(file = here('data', 'contrasts', 'metaD0_contrast1.rds'))
countData <- readRDS(file = here('data', 'contrasts', 'merged_gene_abundance_D0_contrast1.rds'))
