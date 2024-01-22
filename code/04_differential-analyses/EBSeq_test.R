library(EBSeq)
library(here)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# ......................................................................
# LOAD DATA ----
# ......................................................................

meta <- readRDS(file = here('data', 'contrasts', 'metaD0_contrast1.rds'))
#Data file should be genes in rows and samples in columns - Rows with no zeros
countData <- readRDS(file = here('data', 'contrasts', 'merged_gene_abundance_D0_contrast1.rds'))

# ......................................................................
# DIFFERENTIAL ABUNDANCE ----
# ......................................................................

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
countData2 <- countData[rowSums(countData) > 10000,] # ONLY for testing on small batch

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
DE_bact <- annot2[annot2$tax_kingdom == "k__Bacteria",2]
bact_point <- fold[row.names(fold) %in% DE_bact$gene_id,]
DE_metaz <- annot2[annot2$tax_kingdom == "k__Metazoa",2]
metaz_point <- fold[row.names(fold) %in% DE_metaz$gene_id,]
DE_fun <- annot2[grep("mycota", annot2$tax_phylum, ignore.case = TRUE),2]
fun_point <- fold[row.names(fold) %in% DE_fun$gene_id,]
DE_plant <- annot2[annot2$tax_kingdom == "k__Viridiplantae",2]
plant_point <- fold[row.names(fold) %in% DE_plant$gene_id,]
DE_null <- annot2[annot2$tax_kingdom == "k__NULL",2]
null_point <- fold[row.names(fold) %in% DE_null$gene_id,]
DE_other <- annot2[!annot2$gene_id %in% append(append(DE_bact,DE_fun),append(DE_plant, append(DE_null,DE_metaz))),2]
other_point <- fold[row.names(fold) %in% DE_other$gene_id,]


#Plots
#Only plot first 100 000 points for clarity
volcano <- ggplot() +
   geom_point(data = fold[1:100000,], aes(x = mean, y = log2)) +
   geom_point(data = null_point, aes(x = mean, y = log2), color = "purple") +
   geom_point(data = bact_point, aes(x = mean, y = log2), color = "red") +
   geom_point(data = fun_point, aes(x = mean, y = log2), color = "blue") +
   geom_point(data = metaz_point, aes(x = mean, y = log2), color = "yellow") +
   geom_point(data = plant_point, aes(x = mean, y = log2), color = "green") +
   geom_hline(yintercept = 0, linetype = "solid") +
   xlab("Mean read count") + 
   ylab("Log 2 fold change D0 sedmnt OSPW vs CTRL") + 
   scale_x_log10() +
   theme_bw()
volcano 

##### CORRECT FROM HERE ###########

#Sort for exporting
annot2 <- annot2[order(annot2$gene_id),]
DEfound <- data.frame(DE$DEfound)
DEfound <- DEfound[order(DEfound$DE.DEfound),]
fold_DE <- fold[row.names(fold) %in% DE$DEfound,]
fold_DE <- fold_DE[order(row.names(fold_DE)),]
pval_DE <- data.frame(DE$PPMat[row.names(DE$PPMat) %in% DE$DEfound,])
pval_DE <- pval_DE[order(row.names(pval_DE)),]
sum(annot2$gene_id == DEfound) #Sanity check, should be 29 during test
sum(row.names(fold_DE) == DEfound) #Sanity check, should be 29
sum(row.names(pval_DE) == DEfound)#Sanity check, should be 29 
DEfound <- cbind(DEfound,fold_DE, pval_DE, annot2)

#Export
#Matrix of all genes with DE and EE p-values
write.table(as.matrix(DE$PPMat), file = here("output", "PPmat.txt"), eol = "\n", sep = "\t")
#List of the DE transcripts, p-Values, annotations, abundance
write.table(DEfound, file = here("output", "DEfound.txt"), eol = "\n", sep = "\t")
#Status of the transcript: DE, EE, or Filtered+Reason
write.table(as.matrix(DE$Status), file = here("output", "Status.txt"), eol = "\n", sep = "\t")

# ......................................................................
# VARIOUS STATS ----
# ......................................................................


#How many DE gene positive, negative?
length(fun.point.root[fun.point.root$log2 > 0,1]) #23042
length(fun.point.root[fun.point.root$log2 < 0,1]) #232
length(plant.point.root[plant.point.root$log2 > 0,1]) #1295
length(plant.point.root[plant.point.root$log2 < 0,1]) #3061
length(bact_point.root[bact_point.root$log2 > 0,1]) #2231
length(bact_point.root[bact_point.root$log2 < 0,1]) #78
length(other.point.root[other.point.root$log2 > 0,1]) #5895
length(other.point.root[other.point.root$log2 < 0,1]) #863

#How many total transcripts for each taxa ?
length(annot[(annot$tax_kingdom == "Bacteria") & (annot$gene_id %in% row.names(MT.root)),2]) #365 435
length(annot[(annot$tax_kingdom == "Archaea") & (annot$gene_id %in% row.names(MT.root)),2]) #12 792
length(annot[(annot$tax_kingdom == "NULL") & (annot$gene_id %in% row.names(MT.root)),2]) #221 331
length(annot[(annot$tax_kingdom == "Eukaryota") & (annot$gene_id %in% row.names(MT.root)),2]) #533 369
length(annot[(annot$tax_kingdom == "Viruses") & (annot$gene_id %in% row.names(MT.root)),2]) #3 660
length(annot[(annot$tax_phylum == "Streptophyta") & (annot$gene_id %in% row.names(MT.root)),2]) #200 823
length(annot[(annot$tax_phylum == "Ascomycota") & (annot$gene_id %in% row.names(MT.root)),2]) #151 695
length(annot[(annot$tax_phylum == "Basidiomycota") & (annot$gene_id %in% row.names(MT.root)),2]) #46 352
length(annot[(annot$tax_phylum == "Mucoromycota") & (annot$gene_id %in% row.names(MT.root)),2]) #2 186
(533369 - 200823 - 151695 - 46352 - 2186) #132 313

# ......................................................................
# PLOT VARIOUS STATS ----
# ......................................................................

stats.all <- data.frame(
   Taxa = c("Archaea", "Bacteria", "Fungi", "Plants", "Viruses", "unclassified", "Others","Archaea", "Bacteria", "Fungi", "Plants", "Viruses", "unclassified", "Others"),
   Comp = c("Root","Root","Root","Root","Root","Root","Root","Rhizo","Rhizo","Rhizo","Rhizo","Rhizo","Rhizo","Rhizo"),
   Count = c(12792,365435,200233,200823,3660,221331,132313,14943, 430984, 186745,169788, 4255, 254744, 132042)
)
stats.DE <- data.frame(
   Taxa = c("Bacteria-More", "Bacteria-Less","Fungi-More", "Fungi-Less","Plants-More", "Plants-Less","unclassified", "Others","Bacteria-More", "Bacteria-Less","Fungi-More", "Fungi-Less","Plants-More", "Plants-Less","unclassified", "Others"),
   Comp = c("Root","Root","Root","Root","Root","Root","Root","Root", "Rhizo","Rhizo","Rhizo","Rhizo","Rhizo","Rhizo","Rhizo", "Rhizo"),
   Count = c(2231,78,23042,232,1295,3061,5303,6758,7938,6240,149,253,41,178,5224,1583)
)

palette(c("cyan", "red", "blue", "yellow", "green", "purple", "grey"))
stack.stats.all <- ggplot(stats.all, aes(fill = Taxa, y = Count, x = Comp)) + 
   geom_bar( stat = "identity", position = "fill") +
   ylab("Fraction of reads") + 
   scale_fill_manual(values = palette(), guide = guide_legend(label.theme = element_text(face = "italic", size = 8))) +
   theme_bw() +
   scale_y_continuous( expand = c(0,0)) +
   scale_x_discrete(name = "Compartment") +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
stack.stats.all

palette(c("darkred", "red","darkblue", "blue", "yellow","darkgreen", "green", "purple"))
stack.stats.DE <- ggplot(stats.DE, aes(fill = Taxa, y = Count, x = Comp)) + 
   geom_bar( stat = "identity", position = "fill") +
   ylab("Fraction of reads") + 
   scale_fill_manual(values = palette(), guide = guide_legend(label.theme = element_text(face = "italic", size = 8))) +
   theme_bw() +
   scale_y_continuous( expand = c(0,0)) +
   scale_x_discrete(name = "Compartment") +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
stack.stats.DE  


# ......................................................................
# STACK BARCHART TAXA ----
# ......................................................................

```{r stack barchart taxa Rhizo}
#Prokaryotes at phylum
#Create data frame
DEfound.bact <- DEfound[DEfound$tax_kingdom == "Bacteria",]
phylum.table = DEfound.bact %>%
   group_by(tax_phylum) %>%
   summarise(
      n = n(),
      over = sum(log2 > 0, na.rm = TRUE),
      under = sum(log2 < 0, na.rm = TRUE),
      
   ) %>%
   filter(
      n > 100
   )
phylum.table <- phylum.table[-9,] #Remove NULL
row.names(phylum.table) <- phylum.table$tax_phylum
others.phylum <- data.frame("Others", (14178 - sum(phylum.table$n)), (7938 - sum(phylum.table$over)), (6240 - sum(phylum.table$under)))
colnames(others.phylum) <- colnames(phylum.table)
row.names(others.phylum) <- "Others"
phylum.table <- rbind(phylum.table, others.phylum) 
colSums(phylum.table[,2:4]) #Should be 14178

#Get a column for all transcripts
annot.bact <- annot[annot$tax_kingdom == "Bacteria",]
annot.bact <- annot.bact[annot.bact$gene_id %in% row.names(MT),]
MT.bact <- MT[row.names(MT) %in% annot.bact$gene_id,]
annot.bact <- annot.bact[order(annot.bact$gene_id),]
MT.bact <- MT.bact[order(row.names(MT.bact)),]
sum(row.names(MT.bact) == annot.bact$gene_id) #Sanity check
MT.bact.tax <- data.frame(cbind(MT.bact, annot.bact$tax_phylum))

all.column <- MT.bact.tax %>%
   group_by(V13) %>%
   summarise(
      count = n(),
      
      
   )
all.column <- all.column[all.column$V13 %in% phylum.table$tax_phylum,]
row.names(all.column) <- all.column$V13
others.all <- data.frame("Others", (430984 - sum(all.column$count)))
colnames(others.all) <- colnames(all.column)
row.names(others.all) <- "Others"
all.column <- rbind(all.column, others.all) 
sum(all.column$count)

all.column$V13 == phylum.table$tax_phylum #Sanity check
phylum.table <- cbind(phylum.table,all.column$count)
colnames(phylum.table) <- c("Phylum", "n", "More", "Less", "All")

#Prepare for ggplot
phylum.table.long <- gather(phylum.table,subset,counts,3:5) #transform in long format for ggplot
palette(c(brewer.pal(n = 9, name = "Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4", "darkgrey", "white"))
stack.phyla <- ggplot(phylum.table.long, aes(fill = Phylum, y = counts, x = subset)) + 
   geom_bar( stat = "identity", position = "fill") +
   ylab("Fraction of reads") + 
   scale_fill_manual(values = palette(), guide = guide_legend(label.theme = element_text(face = "italic", size = 8))) +
   theme_bw() +
   scale_y_continuous( expand = c(0,0)) +
   scale_x_discrete(name = "Subset") +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
stack.phyla



# ......................................................................
# STACK BARCHART COG ----
# ......................................................................


#COG for FUNGI
#Create data frame
DEfound.root.fun <- DEfound.root[grep(".mycota", DEfound.root$tax_phylum),]
COG.root <- DEfound.root.fun %>%
   group_by(cog_category) %>%
   summarise(
      n = n(),
      over = sum(log2 > 0, na.rm = TRUE),
      under = sum(log2 < 0, na.rm = TRUE),
      
   ) %>%
   filter(
      n > 100
   )
COG.root <- COG.root[-c(4,9,13,14,16,17),] #Remove some categories mixed and NULL... NULL: 14041 13901 140; 
row.names(COG.root) <- COG.root$cog_category
others.cog.root <- data.frame("Others", (23274 - 14041 - sum(COG.root$n)), (23042 - 13901 - sum(COG.root$over)), (232 - 140 - um(COG.root$under))) #Create other bu omit NULL
colnames(others.cog.root) <- colnames(COG.root)
row.names(others.cog.root) = "Others"
COG.root <- rbind(COG.root, others.cog.root) 
colSums(COG.root[,2:4])

#Get a column for all fungal transcripts
annot.fun <- annot[grep(".mycota", annot$tax_phylum),]
annot.fun.root <- annot.fun[annot.fun$gene_id %in% row.names(MT.root),]
MT.root.fun <- MT.root[row.names(MT.root) %in% annot.fun$gene_id,]
annot.fun.root <- annot.fun.root[order(annot.fun.root$gene_id),]
MT.root.fun <- MT.root.fun[order(row.names(MT.root.fun)),]
sum(row.names(MT.root.fun) == annot.fun.root$gene_id)
MT.root.fun.cog <- data.frame(cbind(MT.root.fun, annot.fun.root$cog_category))

COG.root.all = MT.root.fun.cog %>%
   group_by(V12) %>%
   summarise(
      count = n(),
      
      
   )
COG.root.all <- COG.root.all[COG.root.all$V12 %in% COG.root$cog_category,]
row.names(COG.root.all) <- COG.root.all$V12
others.cog.root.all <- data.frame("Others", (200233 - 134958 - sum(COG.root.all$count))) #Omit NULL
colnames(others.cog.root.all) <- colnames(COG.root.all)
row.names(others.cog.root.all) <- "Others"
COG.root.all <- rbind(COG.root.all, others.cog.root.all) 
sum(COG.root.all$count)

COG.root.all$V12 == COG.root$cog_category
COG.root <- cbind(COG.root,COG.root.all$count)
colnames(COG.root) <- c("COG_category", "n", "More", "Less", "All")

#Prepare for ggplot
COG.root.long <- gather(COG.root,subset,counts,3:5) #transform in long format for ggplot
palette(c(brewer.pal(n = 9, name = "Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4", "brown3", "cyan", "pink"))
stack.COG.root <- ggplot(COG.root.long, aes(fill = COG_category, y = counts, x = subset)) + 
   geom_bar( stat = "identity", position = "fill") +
   ylab("Fraction of reads") + 
   scale_fill_manual(values = palette(), guide = guide_legend(label.theme = element_text(face = "plain", size = 8))) +
   theme_bw() +
   scale_y_continuous( expand = c(0,0)) +
   scale_x_discrete(name = "Subset") +
   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
stack.COG.root


# ......................................................................
# TOP 50 DE GENES  ----
# ......................................................................


###ROOTS
#Sort based on P-value of EE
DEfound.root.50 <- DEfound.root[DEfound.root$PPEE == 0,] #lots of Pvalue = 0... 1491. Then take most abundant.
DEfound.root.50 <- DEfound.root.50[order(DEfound.root.50$mean, decreasing = TRUE),]
DEfound.root.50.defined <- DEfound.root.50[DEfound.root.50$cog_function != "NULL",] #Remove undefined
DEfound.root.50.defined <- DEfound.root.50.defined[1:50, c(1:4, 24,25,34,37)] #Keep omly relevant columns
#Export
write.table(DEfound.root.50.defined, file = here("output", "Table1.txt"), eol = "\n", sep = "\t")

###RHIZO
#Sort based on P-value of EE
DEfound.50 <- DEfound[order(DEfound$PPEE),]
DEfound.50.defined <- DEfound.50[DEfound.50$cog_function != "NULL",] #Remove undefined
DEfound.50.defined <- DEfound.50.defined[1:100, c(1:4, 24,25,34,37)] #Keep only relevant columns, keep 100, because lot of irrelevant orgs; need to manually remove
#Export
write.table(DEfound.50.defined, file = here("output", "Table2.txt"), eol = "\n", sep = "\t")



# ......................................................................
# ANOVA TAXA ----
# ......................................................................

```{r anova taxa Rhizo}
#Prokaryotes at phylum
#Get table with all reps
MT.bact.tax.found <- MT.bact.tax[row.names(MT.bact.tax) %in% row.names(DEfound.bact),] #Get only DEfound
DEfound.bact.s <- DEfound.bact[order(row.names(DEfound.bact)),] #Sort
sum(row.names(DEfound.bact.s) == row.names(MT.bact.tax.found)) #14178
MT.bact.tax.found <- data.frame(cbind(MT.bact.tax.found, DEfound.bact.s$log2))#Add column for log2
MT.bact.tax.found.over <- MT.bact.tax.found[DEfound.bact.s$log2 > 0,]
MT.bact.tax.found.under <- MT.bact.tax.found[DEfound.bact.s$log2 < 0,]
MT.bact.tax.found.over <- data.frame(apply(MT.bact.tax.found.over[,1:12], 2,      # Specify own function within apply
                                                 function(x) as.numeric(as.character(x))), MT.bact.tax.found.over[,13:14])
MT.bact.tax.found.under <- data.frame(apply(MT.bact.tax.found.under[,1:12], 2,      # Specify own function within apply
                                                  function(x) as.numeric(as.character(x))), MT.bact.tax.found.under[,13:14])
MT.bact.tax <- data.frame(apply(MT.bact.tax[,1:12], 2,      # Specify own function within apply
                                      function(x) as.numeric(as.character(x))), MT.bact.tax[,13])



#DEFound, only group_by phylum, keep all samples
summary.bact.over <- MT.bact.tax.found.over %>%
   group_by(V13) %>%
   summarise(
      B1P2S.over = sum(B1P2S),
      B1P3S.over = sum(B1P3S),
      B2P1S.over = sum(B2P1S),
      B2P2S.over = sum(B2P2S),
      B3P2S.over = sum(B3P2S),
      B3P8S.over = sum(B3P8S),
      B4P5S.over = sum(B4P5S),
      B4P7S.over = sum(B4P7S),
      B5P2S.over = sum(B5P2S),
      B5P3S.over = sum(B5P3S),
      B6P2S.over = sum(B6P2S),
      B6P6S.over = sum(B6P6S),
      
   )
summary.bact.under <- MT.bact.tax.found.under %>%
   group_by(V13) %>%
   summarise(
      B1P2S.under = sum(B1P2S),
      B1P3S.under = sum(B1P3S),
      B2P1S.under = sum(B2P1S),
      B2P2S.under = sum(B2P2S),
      B3P2S.under = sum(B3P2S),
      B3P8S.under = sum(B3P8S),
      B4P5S.under = sum(B4P5S),
      B4P7S.under = sum(B4P7S),
      B5P2S.under = sum(B5P2S),
      B5P3S.under = sum(B5P3S),
      B6P2S.under = sum(B6P2S),
      B6P6S.under = sum(B6P6S),
      
   )

summary.bact.all <- MT.bact.tax %>%
   group_by(MT.bact.tax...13.) %>%
   summarise(
      B1P2S.all = sum(B1P2S),
      B1P3S.all = sum(B1P3S),
      B2P1S.all = sum(B2P1S),
      B2P2S.all = sum(B2P2S),
      B3P2S.all = sum(B3P2S),
      B3P8S.all = sum(B3P8S),
      B4P5S.all = sum(B4P5S),
      B4P7S.all = sum(B4P7S),
      B5P2S.all = sum(B5P2S),
      B5P3S.all = sum(B5P3S),
      B6P2S.all = sum(B6P2S),
      B6P6S.all = sum(B6P6S),
      
   )

#Put everything in relative
#All
summary.bact.all$Phylum <- summary.bact.all$MT.bact.tax...13.
summary.bact.all <- summary.bact.all[,-1]
summary.bact.all.rel <- data.frame(t(apply(summary.bact.all[,1:12], 1, "/", colSums(summary.bact.all[,1:12]))))
colSums(summary.bact.all.rel)
summary.bact.all.rel$Phylum <- summary.bact.all$Phylum
#Over
summary.bact.over$Phylum <- summary.bact.over$V13
summary.bact.over <- summary.bact.over[,-1]
summary.bact.over.rel <- data.frame(t(apply(summary.bact.over[,1:12], 1, "/", colSums(summary.bact.over[,1:12]))))
colSums(summary.bact.over.rel)
summary.bact.over.rel$Phylum <- summary.bact.over$Phylum
#under
summary.bact.under$Phylum <- summary.bact.under$V13
summary.bact.under <- summary.bact.under[,-1]
summary.bact.under.rel <- data.frame(t(apply(summary.bact.under[,1:12], 1, "/", colSums(summary.bact.under[,1:12]))))
colSums(summary.bact.under.rel)
summary.bact.under.rel$Phylum <- summary.bact.under$Phylum

#ALL vs. OVER
#Acidobacteria
"All:"; summary.bact.all.rel[2,13]; mean(as.numeric(summary.bact.all.rel[2,1:12]))
"Over:"; summary.bact.over.rel[2,13]; mean(as.numeric(summary.bact.over.rel[2,1:12]))
t.test(t(summary.bact.all.rel[2,1:12]), t(summary.bact.over.rel[2,1:12]), paired = TRUE)
#Actinobacteria
"All:"; summary.bact.all.rel[3,13]; mean(as.numeric(summary.bact.all.rel[3,1:12]))
"Over:"; summary.bact.over.rel[3,13]; mean(as.numeric(summary.bact.over.rel[3,1:12]))
t.test(t(summary.bact.all.rel[3,1:12]), t(summary.bact.over.rel[3,1:12]), paired = TRUE)
#Bacteroidetes
"All:"; summary.bact.all.rel[6,13]; mean(as.numeric(summary.bact.all.rel[6,1:12]))
"Over:"; summary.bact.over.rel[5,13]; mean(as.numeric(summary.bact.over.rel[5,1:12]))
t.test(t(summary.bact.all.rel[6,1:12]), t(summary.bact.over.rel[5,1:12]), paired = TRUE)
#Candidatus Rokubacteria
"All:"; summary.bact.all.rel[91,13]; mean(as.numeric(summary.bact.all.rel[91,1:12]))
"Over:"; summary.bact.over.rel[20,13]; mean(as.numeric(summary.bact.over.rel[20,1:12]))
t.test(t(summary.bact.all.rel[91,1:12]), t(summary.bact.over.rel[20,1:12]), paired = TRUE)
#Chloroflexi
"All:"; summary.bact.all.rel[115,13]; mean(as.numeric(summary.bact.all.rel[115,1:12]))
"Over:"; summary.bact.over.rel[24,13]; mean(as.numeric(summary.bact.over.rel[24,1:12]))
t.test(t(summary.bact.all.rel[115,1:12]), t(summary.bact.over.rel[24,1:12]), paired = TRUE)
#Cyanobacteria
"All:"; summary.bact.all.rel[118,13]; mean(as.numeric(summary.bact.all.rel[118,1:12]))
"Over:"; summary.bact.over.rel[25,13]; mean(as.numeric(summary.bact.over.rel[25,1:12]))
t.test(t(summary.bact.all.rel[118,1:12]), t(summary.bact.over.rel[25,1:12]), paired = TRUE)
#Gemmatimonadetes
"All:"; summary.bact.all.rel[126,13]; mean(as.numeric(summary.bact.all.rel[126,1:12]))
"Over:"; summary.bact.over.rel[29,13]; mean(as.numeric(summary.bact.over.rel[29,1:12]))
t.test(t(summary.bact.all.rel[126,1:12]), t(summary.bact.over.rel[29,1:12]), paired = TRUE)
#Nitrospira
"All:"; summary.bact.all.rel[131,13]; mean(as.numeric(summary.bact.all.rel[131,1:12]))
"Over:"; summary.bact.over.rel[33,13]; mean(as.numeric(summary.bact.over.rel[33,1:12]))
t.test(t(summary.bact.all.rel[131,1:12]), t(summary.bact.over.rel[33,1:12]), paired = TRUE)
#Planctomycetes
"All:"; summary.bact.all.rel[133,13]; mean(as.numeric(summary.bact.all.rel[133,1:12]))
"Over:"; summary.bact.over.rel[35,13]; mean(as.numeric(summary.bact.over.rel[35,1:12]))
t.test(t(summary.bact.all.rel[133,1:12]), t(summary.bact.over.rel[35,1:12]), paired = TRUE)
#Proteobacteria
"All:"; summary.bact.all.rel[134,13]; mean(as.numeric(summary.bact.all.rel[134,1:12]))
"Over:"; summary.bact.over.rel[36,13]; mean(as.numeric(summary.bact.over.rel[36,1:12]))
t.test(t(summary.bact.all.rel[134,1:12]), t(summary.bact.over.rel[36,1:12]), paired = TRUE)
#Verrucomicrobia
"All:"; summary.bact.all.rel[141,13]; mean(as.numeric(summary.bact.all.rel[141,1:12]))
"Over:"; summary.bact.over.rel[39,13]; mean(as.numeric(summary.bact.over.rel[39,1:12]))
t.test(t(summary.bact.all.rel[141,1:12]), t(summary.bact.over.rel[39,1:12]), paired = TRUE)

#ALL vs. UNDER
#Acidobacteria
"All:"; summary.bact.all.rel[2,13]; mean(as.numeric(summary.bact.all.rel[2,1:12]))
"under:"; summary.bact.under.rel[2,13]; mean(as.numeric(summary.bact.under.rel[2,1:12]))
t.test(t(summary.bact.all.rel[2,1:12]), t(summary.bact.under.rel[2,1:12]), paired = TRUE)
#Actinobacteria
"All:"; summary.bact.all.rel[3,13]; mean(as.numeric(summary.bact.all.rel[3,1:12]))
"under:"; summary.bact.under.rel[3,13]; mean(as.numeric(summary.bact.under.rel[3,1:12]))
t.test(t(summary.bact.all.rel[3,1:12]), t(summary.bact.under.rel[3,1:12]), paired = TRUE)
#Bacteroidetes
"All:"; summary.bact.all.rel[6,13]; mean(as.numeric(summary.bact.all.rel[6,1:12]))
"under:"; summary.bact.under.rel[6,13]; mean(as.numeric(summary.bact.under.rel[6,1:12]))
t.test(t(summary.bact.all.rel[6,1:12]), t(summary.bact.under.rel[6,1:12]), paired = TRUE)
#Candidatus Rokubacteria
"All:"; summary.bact.all.rel[91,13]; mean(as.numeric(summary.bact.all.rel[91,1:12]))
"under:"; summary.bact.under.rel[35,13]; mean(as.numeric(summary.bact.under.rel[35,1:12]))
t.test(t(summary.bact.all.rel[91,1:12]), t(summary.bact.under.rel[35,1:12]), paired = TRUE)
#Chloroflexi
"All:"; summary.bact.all.rel[115,13]; mean(as.numeric(summary.bact.all.rel[115,1:12]))
"under:"; summary.bact.under.rel[46,13]; mean(as.numeric(summary.bact.under.rel[46,1:12]))
t.test(t(summary.bact.all.rel[115,1:12]), t(summary.bact.under.rel[46,1:12]), paired = TRUE)
#Cyanobacteria
"All:"; summary.bact.all.rel[118,13]; mean(as.numeric(summary.bact.all.rel[118,1:12]))
"under:"; summary.bact.under.rel[47,13]; mean(as.numeric(summary.bact.under.rel[47,1:12]))
t.test(t(summary.bact.all.rel[118,1:12]), t(summary.bact.under.rel[47,1:12]), paired = TRUE)
#Gemmatimonadetes
"All:"; summary.bact.all.rel[126,13]; mean(as.numeric(summary.bact.all.rel[126,1:12]))
"under:"; summary.bact.under.rel[52,13]; mean(as.numeric(summary.bact.under.rel[52,1:12]))
t.test(t(summary.bact.all.rel[126,1:12]), t(summary.bact.under.rel[52,1:12]), paired = TRUE)
#Nitrospira
"All:"; summary.bact.all.rel[131,13]; mean(as.numeric(summary.bact.all.rel[131,1:12]))
"under:"; summary.bact.under.rel[56,13]; mean(as.numeric(summary.bact.under.rel[56,1:12]))
t.test(t(summary.bact.all.rel[131,1:12]), t(summary.bact.under.rel[56,1:12]), paired = TRUE)
#Planctomycetes
"All:"; summary.bact.all.rel[133,13]; mean(as.numeric(summary.bact.all.rel[133,1:12]))
"under:"; summary.bact.under.rel[58,13]; mean(as.numeric(summary.bact.under.rel[58,1:12]))
t.test(t(summary.bact.all.rel[133,1:12]), t(summary.bact.under.rel[58,1:12]), paired = TRUE)
#Proteobacteria
"All:"; summary.bact.all.rel[134,13]; mean(as.numeric(summary.bact.all.rel[134,1:12]))
"under:"; summary.bact.under.rel[59,13]; mean(as.numeric(summary.bact.under.rel[59,1:12]))
t.test(t(summary.bact.all.rel[134,1:12]), t(summary.bact.under.rel[59,1:12]), paired = TRUE)
#Verrucomicrobia
"All:"; summary.bact.all.rel[141,13]; mean(as.numeric(summary.bact.all.rel[141,1:12]))
"under:"; summary.bact.under.rel[65,13]; mean(as.numeric(summary.bact.under.rel[65,1:12]))
t.test(t(summary.bact.all.rel[141,1:12]), t(summary.bact.under.rel[65,1:12]), paired = TRUE)





# ......................................................................
# ANOVA COG  ----
# ......................................................................


#Prokaryotes at phylum
#Get table with all reps
MT.bact.cog.found <- MT.bact.cog[row.names(MT.bact.cog) %in% row.names(DEfound.bact),] #Get only DEfound
DEfound.bact.s <- DEfound.bact[order(row.names(DEfound.bact)),] #Sort
sum(row.names(DEfound.bact.s) == row.names(MT.bact.cog.found)) #14178
MT.bact.cog.found <- data.frame(cbind(MT.bact.cog.found, DEfound.bact.s$log2))#Add column for log2
MT.bact.cog.found.over <- MT.bact.cog.found[DEfound.bact.s$log2 > 0,]
MT.bact.cog.found.under <- MT.bact.cog.found[DEfound.bact.s$log2 < 0,]
MT.bact.cog.found.over <- data.frame(apply(MT.bact.cog.found.over[,1:12], 2,  # Specify own function within apply
                                                 function(x) as.numeric(as.character(x))), MT.bact.cog.found.over[,13:14])
MT.bact.cog.found.under <- data.frame(apply(MT.bact.cog.found.under[,1:12], 2,  # Specify own function within apply
                                                  function(x) as.numeric(as.character(x))), MT.bact.cog.found.under[,13:14])
MT.bact.cog <- data.frame(apply(MT.bact.cog[,1:12], 2, # Specify own function within apply
                                      function(x) as.numeric(as.character(x))), MT.bact.cog[,13])



#DEFound, only group_by phylum, keep all samples
summary.cog.over <- MT.bact.cog.found.over %>%
   group_by(V13) %>%
   summarise(
      B1P2S.over = sum(B1P2S),
      B1P3S.over = sum(B1P3S),
      B2P1S.over = sum(B2P1S),
      B2P2S.over = sum(B2P2S),
      B3P2S.over = sum(B3P2S),
      B3P8S.over = sum(B3P8S),
      B4P5S.over = sum(B4P5S),
      B4P7S.over = sum(B4P7S),
      B5P2S.over = sum(B5P2S),
      B5P3S.over = sum(B5P3S),
      B6P2S.over = sum(B6P2S),
      B6P6S.over = sum(B6P6S),
      
   )
summary.cog.under <- MT.bact.cog.found.under %>%
   group_by(V13) %>%
   summarise(
      B1P2S.under = sum(B1P2S),
      B1P3S.under = sum(B1P3S),
      B2P1S.under = sum(B2P1S),
      B2P2S.under = sum(B2P2S),
      B3P2S.under = sum(B3P2S),
      B3P8S.under = sum(B3P8S),
      B4P5S.under = sum(B4P5S),
      B4P7S.under = sum(B4P7S),
      B5P2S.under = sum(B5P2S),
      B5P3S.under = sum(B5P3S),
      B6P2S.under = sum(B6P2S),
      B6P6S.under = sum(B6P6S),
      
   )

summary.cog.all <- MT.bact.cog %>%
   group_by(MT.bact.cog...13.) %>%
   summarise(
      B1P2S.all = sum(B1P2S),
      B1P3S.all = sum(B1P3S),
      B2P1S.all = sum(B2P1S),
      B2P2S.all = sum(B2P2S),
      B3P2S.all = sum(B3P2S),
      B3P8S.all = sum(B3P8S),
      B4P5S.all = sum(B4P5S),
      B4P7S.all = sum(B4P7S),
      B5P2S.all = sum(B5P2S),
      B5P3S.all = sum(B5P3S),
      B6P2S.all = sum(B6P2S),
      B6P6S.all = sum(B6P6S),
      
   )

#Put everything in relative
#All
summary.cog.all$COG_cat <- summary.cog.all$MT.bact.cog...13.
summary.cog.all <- summary.cog.all[,-1]
summary.cog.all.rel <- data.frame(t(apply(summary.cog.all[,1:12], 1, "/", colSums(summary.cog.all[,1:12]))))
colSums(summary.cog.all.rel)
summary.cog.all.rel$COG_cat <- summary.cog.all$COG_cat
#Over
summary.cog.over$COG_cat <- summary.cog.over$V13
summary.cog.over <- summary.cog.over[,-1]
summary.cog.over.rel <- data.frame(t(apply(summary.cog.over[,1:12], 1, "/", colSums(summary.cog.over[,1:12]))))
colSums(summary.cog.over.rel)
summary.cog.over.rel$COG_cat <- summary.cog.over$COG_cat
#under
summary.cog.under$COG_cat <- summary.cog.under$V13
summary.cog.under <- summary.cog.under[,-1]
summary.cog.under.rel <- data.frame(t(apply(summary.cog.under[,1:12], 1, "/", colSums(summary.cog.under[,1:12]))))
colSums(summary.cog.under.rel)
summary.cog.under.rel$COG_cat <- summary.cog.under$COG_cat

#ALL vs. OVER
#AA trans and metab
"All:"; summary.cog.all.rel[6,13]; mean(as.numeric(summary.cog.all.rel[6,1:12]))
"Over:"; summary.cog.over.rel[3,13]; mean(as.numeric(summary.cog.over.rel[3,1:12]))
t.test(t(summary.cog.all.rel[6,1:12]), t(summary.cog.over.rel[3,1:12]), paired = TRUE)
#Carb trans and metabolism
"All:"; summary.cog.all.rel[17,13]; mean(as.numeric(summary.cog.all.rel[17,1:12]))
"Over:"; summary.cog.over.rel[9,13]; mean(as.numeric(summary.cog.over.rel[9,1:12]))
t.test(t(summary.cog.all.rel[17,1:12]), t(summary.cog.over.rel[9,1:12]), paired = TRUE)
#Cell envelope
"All:"; summary.cog.all.rel[28,13]; mean(as.numeric(summary.cog.all.rel[28,1:12]))
"Over:"; summary.cog.over.rel[16,13]; mean(as.numeric(summary.cog.over.rel[16,1:12]))
t.test(t(summary.cog.all.rel[28,1:12]), t(summary.cog.over.rel[16,1:12]), paired = TRUE)
#Cell motility and secretion
"All:"; summary.cog.all.rel[36,13]; mean(as.numeric(summary.cog.all.rel[36,1:12]))
"Over:"; summary.cog.over.rel[21,13]; mean(as.numeric(summary.cog.over.rel[21,1:12]))
t.test(t(summary.cog.all.rel[36,1:12]), t(summary.cog.over.rel[21,1:12]), paired = TRUE)
#Coenzyme metabolism
"All:"; summary.cog.all.rel[45,13]; mean(as.numeric(summary.cog.all.rel[45,1:12]))
"Over:"; summary.cog.over.rel[25,13]; mean(as.numeric(summary.cog.over.rel[25,1:12]))
t.test(t(summary.cog.all.rel[45,1:12]), t(summary.cog.over.rel[25,1:12]), paired = TRUE)
#DNA replication
"All:"; summary.cog.all.rel[58,13]; mean(as.numeric(summary.cog.all.rel[58,1:12]))
"Over:"; summary.cog.over.rel[32,13]; mean(as.numeric(summary.cog.over.rel[32,1:12]))
t.test(t(summary.cog.all.rel[58,1:12]), t(summary.cog.over.rel[32,1:12]), paired = TRUE)
#Energy production and conversion
"All:"; summary.cog.all.rel[64,13]; mean(as.numeric(summary.cog.all.rel[64,1:12]))
"Over:"; summary.cog.over.rel[35,13]; mean(as.numeric(summary.cog.over.rel[35,1:12]))
t.test(t(summary.cog.all.rel[64,1:12]), t(summary.cog.over.rel[35,1:12]), paired = TRUE)
#Inorganic ion
"All:"; summary.cog.all.rel[77,13]; mean(as.numeric(summary.cog.all.rel[77,1:12]))
"Over:"; summary.cog.over.rel[42,13]; mean(as.numeric(summary.cog.over.rel[42,1:12]))
t.test(t(summary.cog.all.rel[77,1:12]), t(summary.cog.over.rel[42,1:12]), paired = TRUE)
#Intracellular trafficking
"All:"; summary.cog.all.rel[83,13]; mean(as.numeric(summary.cog.all.rel[83,1:12]))
"Over:"; summary.cog.over.rel[45,13]; mean(as.numeric(summary.cog.over.rel[45,1:12]))
t.test(t(summary.cog.all.rel[83,1:12]), t(summary.cog.over.rel[45,1:12]), paired = TRUE)
#Lipid
"All:"; summary.cog.all.rel[89,13]; mean(as.numeric(summary.cog.all.rel[89,1:12]))
"Over:"; summary.cog.over.rel[48,13]; mean(as.numeric(summary.cog.over.rel[48,1:12]))
t.test(t(summary.cog.all.rel[89,1:12]), t(summary.cog.over.rel[48,1:12]), paired = TRUE)
#Posttranslational
"All:"; summary.cog.all.rel[102,13]; mean(as.numeric(summary.cog.all.rel[102,1:12]))
"Over:"; summary.cog.over.rel[55,13]; mean(as.numeric(summary.cog.over.rel[55,1:12]))
t.test(t(summary.cog.all.rel[102,1:12]), t(summary.cog.over.rel[55,1:12]), paired = TRUE)
#Signal
"All:"; summary.cog.all.rel[125,13]; mean(as.numeric(summary.cog.all.rel[125,1:12]))
"Over:"; summary.cog.over.rel[64,13]; mean(as.numeric(summary.cog.over.rel[64,1:12]))
t.test(t(summary.cog.all.rel[125,1:12]), t(summary.cog.over.rel[64,1:12]), paired = TRUE)
#Transcription
"All:"; summary.cog.all.rel[129,13]; mean(as.numeric(summary.cog.all.rel[129,1:12]))
"Over:"; summary.cog.over.rel[67,13]; mean(as.numeric(summary.cog.over.rel[67,1:12]))
t.test(t(summary.cog.all.rel[129,1:12]), t(summary.cog.over.rel[67,1:12]), paired = TRUE)
#Translation
"All:"; summary.cog.all.rel[135,13]; mean(as.numeric(summary.cog.all.rel[135,1:12]))
"Over:"; summary.cog.over.rel[72,13]; mean(as.numeric(summary.cog.over.rel[72,1:12]))
t.test(t(summary.cog.all.rel[135,1:12]), t(summary.cog.over.rel[72,1:12]), paired = TRUE)

#ALL vs. UNDER
#AA trans and metab
"All:"; summary.cog.all.rel[6,13]; mean(as.numeric(summary.cog.all.rel[6,1:12]))
"Under:"; summary.cog.under.rel[4,13]; mean(as.numeric(summary.cog.under.rel[4,1:12]))
t.test(t(summary.cog.all.rel[6,1:12]), t(summary.cog.under.rel[4,1:12]), paired = TRUE)
#Carb trans and metabolism
"All:"; summary.cog.all.rel[17,13]; mean(as.numeric(summary.cog.all.rel[17,1:12]))
"Under:"; summary.cog.under.rel[13,13]; mean(as.numeric(summary.cog.under.rel[13,1:12]))
t.test(t(summary.cog.all.rel[17,1:12]), t(summary.cog.under.rel[13,1:12]), paired = TRUE)
#Cell envelope
"All:"; summary.cog.all.rel[28,13]; mean(as.numeric(summary.cog.all.rel[28,1:12]))
"Under:"; summary.cog.under.rel[16,13]; mean(as.numeric(summary.cog.under.rel[16,1:12]))
t.test(t(summary.cog.all.rel[28,1:12]), t(summary.cog.under.rel[16,1:12]), paired = TRUE)
#Cell motility and secretion
"All:"; summary.cog.all.rel[36,13]; mean(as.numeric(summary.cog.all.rel[36,1:12]))
"Under:"; summary.cog.under.rel[20,13]; mean(as.numeric(summary.cog.under.rel[20,1:12]))
t.test(t(summary.cog.all.rel[36,1:12]), t(summary.cog.under.rel[20,1:12]), paired = TRUE)
#Coenzyme metabolism
"All:"; summary.cog.all.rel[45,13]; mean(as.numeric(summary.cog.all.rel[45,1:12]))
"Under:"; summary.cog.under.rel[25,13]; mean(as.numeric(summary.cog.under.rel[25,1:12]))
t.test(t(summary.cog.all.rel[45,1:12]), t(summary.cog.under.rel[25,1:12]), paired = TRUE)
#DNA replication
"All:"; summary.cog.all.rel[58,13]; mean(as.numeric(summary.cog.all.rel[58,1:12]))
"Under:"; summary.cog.under.rel[30,13]; mean(as.numeric(summary.cog.under.rel[30,1:12]))
t.test(t(summary.cog.all.rel[58,1:12]), t(summary.cog.under.rel[30,1:12]), paired = TRUE)
#Energy production and conversion
"All:"; summary.cog.all.rel[64,13]; mean(as.numeric(summary.cog.all.rel[64,1:12]))
"Under:"; summary.cog.under.rel[35,13]; mean(as.numeric(summary.cog.under.rel[35,1:12]))
t.test(t(summary.cog.all.rel[64,1:12]), t(summary.cog.under.rel[35,1:12]), paired = TRUE)
#Inorganic ion
"All:"; summary.cog.all.rel[77,13]; mean(as.numeric(summary.cog.all.rel[77,1:12]))
"Under:"; summary.cog.under.rel[41,13]; mean(as.numeric(summary.cog.under.rel[41,1:12]))
t.test(t(summary.cog.all.rel[77,1:12]), t(summary.cog.under.rel[41,1:12]), paired = TRUE)
#Intracellular trafficking
"All:"; summary.cog.all.rel[83,13]; mean(as.numeric(summary.cog.all.rel[83,1:12]))
"Under:"; summary.cog.under.rel[44,13]; mean(as.numeric(summary.cog.under.rel[44,1:12]))
t.test(t(summary.cog.all.rel[83,1:12]), t(summary.cog.under.rel[44,1:12]), paired = TRUE)
#Lipid
"All:"; summary.cog.all.rel[89,13]; mean(as.numeric(summary.cog.all.rel[89,1:12]))
"Under:"; summary.cog.under.rel[47,13]; mean(as.numeric(summary.cog.under.rel[47,1:12]))
t.test(t(summary.cog.all.rel[89,1:12]), t(summary.cog.under.rel[47,1:12]), paired = TRUE)
#Posttranslational
"All:"; summary.cog.all.rel[102,13]; mean(as.numeric(summary.cog.all.rel[102,1:12]))
"Under:"; summary.cog.under.rel[53,13]; mean(as.numeric(summary.cog.under.rel[53,1:12]))
t.test(t(summary.cog.all.rel[102,1:12]), t(summary.cog.under.rel[53,1:12]), paired = TRUE)
#Signal
"All:"; summary.cog.all.rel[125,13]; mean(as.numeric(summary.cog.all.rel[125,1:12]))
"Under:"; summary.cog.under.rel[63,13]; mean(as.numeric(summary.cog.under.rel[63,1:12]))
t.test(t(summary.cog.all.rel[125,1:12]), t(summary.cog.under.rel[63,1:12]), paired = TRUE)
#Transcription
"All:"; summary.cog.all.rel[129,13]; mean(as.numeric(summary.cog.all.rel[129,1:12]))
"Under:"; summary.cog.under.rel[65,13]; mean(as.numeric(summary.cog.under.rel[65,1:12]))
t.test(t(summary.cog.all.rel[129,1:12]), t(summary.cog.under.rel[65,1:12]), paired = TRUE)
#Translation
"All:"; summary.cog.all.rel[135,13]; mean(as.numeric(summary.cog.all.rel[135,1:12]))
"Under:"; summary.cog.under.rel[68,13]; mean(as.numeric(summary.cog.under.rel[68,1:12]))
t.test(t(summary.cog.all.rel[135,1:12]), t(summary.cog.under.rel[68,1:12]), paired = TRUE)


# ......................................................................
# ANOVA COG  ----
# ......................................................................


