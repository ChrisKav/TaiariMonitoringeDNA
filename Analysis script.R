library(readxl)
library(openxlsx)
library(splitstackshape)
library(data.table)
library(lubridate)
library(tidyverse)
library(taxa)
library(phyloseq)
library(metacoder)
library(dplyr)

data <- read_excel("eDNA survey 2022.xlsx")
coa <- read_xlsx("WLJ602811.xlsx", sheet=1)
coa <- coa[19:nrow(coa),]
colnames(coa) <- coa[1,]
coa <- data.frame(coa[2:nrow(coa),])
taxa_data <- read_xlsx("WLJ602811.xlsx", sheet=2)
seq_data <- read_xlsx("WLJ602811.xlsx", sheet=3)
TICI <- read_xlsx("WLJ602811.xlsx", sheet=4)

df <- cSplit(data, "UID", sep = ",", direction = "long")

df$`Retrieval date` <- as.Date(df$`Retrieval date`, format = "%y/%m/%d")
df$`Deployment date` <- as.Date(df$`Deployment date`, format = "%y/%m/%d")

Retrieval = data.frame(
  date=df$`Retrieval date`,
  time=format(df$`Retrieval time`, "%H:%M")
)

Deployment = data.frame(
  date=df$`Deployment date`,
  time=format(df$`Deployment time`, "%H:%M")
)

Deployment <- as.POSIXct(paste(Deployment$date, Deployment$time), format="%Y-%m-%d %H:%M")
Retrieval <- as.POSIXct(paste(Retrieval$date, Retrieval$time), format="%Y-%m-%d %H:%M")

df$period <- round(difftime(Retrieval, Deployment, units = "hours"),1)

df <- df[which(df$UID %in% coa$UID),]

write.xlsx(df, "COA list.xlsx")

### OTU Table
otumat <- data.frame(taxa_data[c(3,6:193)])

#Need to change TaxID to match ncbi database CHECK!!!
otumat[which(otumat$TaxID == 10000005),1] = 27462 # Austrosimulium changed to genus
otumat[which(otumat$TaxID == 10000018),1] = 126351 # Changed to Gal spD
otumat[which(otumat$TaxID == 10000020),1] = 89553 # Changed to Giant Kokopu
otumat[which(otumat$TaxID == 10000033),1] = 309669 # Changed to genus
otumat[which(otumat$TaxID == 10000037),1] = 65076 # Changed to genus
otumat[which(otumat$TaxID == 10000038),1] = 226931 # Changed to Common Bully
otumat[which(otumat$TaxID == 10000043),1] = 98305 #Changed to genus
otumat[which(otumat$TaxID == 10000052),1] = 75837 # Changed to Brown Teal
otumat[which(otumat$TaxID == 10000053),1] = 126303 # Changed to Gal dep
otumat[which(otumat$TaxID == 10000054),1] = 126303 # Changed to Gal dep

#Combine rows with same TaxID
otumat <- aggregate(x = otumat, by = list(otumat$TaxID), FUN = function(x) na.omit(x)[1])[,-1]

#write.csv(otumat, "x1.csv") # for checking

row.names(otumat) <- otumat$TaxID
otumat$TaxID <- NULL
otumat <- as.matrix(otumat)
OTU = otu_table(otumat, taxa_are_rows = TRUE)

### USE library(txa and metacoder) to do taxonomy for phyloseq
# x <- lookup_tax_data(taxdf$ID, type = "taxon_id") #takes a while to run
y <- taxonomy_table(x, use_ranks = c("kingdom", "phylum", "class", "order",
        "family", "genus", "species"), add_id_col = TRUE)
y <- data.frame(y)
rownames(y) <- y$taxon_id
y$taxon_id <- NULL
taxmat <- as.matrix(y)

# write.csv(y, "tax_table.csv") # Check

TAX = tax_table(taxmat)

### sample_datae construction
sampledata <- data.frame(cbind(df$UID, as.Date(df$`Deployment date`), df$period, df$`Waterway name`, df$Site, TICI$TICI_value, TICI$TICI_rating, TICI$TICI_reliability))
row.names(sampledata) <- df$UID
colnames(sampledata) <- c('X.SampleID', 'Date', 'Period', 'Waterway', 'Site', 'TICI_val', 'TICI_rating', 'TICI_rel')
sampledata <- sample_data(sampledata)
sample_names(sampledata) <- paste0("X", sample_names(sampledata))

#Build phyloseq object
physeq = phyloseq(OTU, TAX, sampledata) # Sample data not matching UID
physeq

#can use metacoder to revert phyloseq object and produce taxmaps for each segment of river
#Need to get date in correct format

### Merge samples to sample site

phymer = merge_samples(physeq, "Site")

#need to build usable sample_data()

#ordination plot of all site
all_pcoa <- ordinate(
  physeq = phymer, 
  method = "PCoA", 
  distance = "bray"
)

plot_ordination(
  physeq = phymer,                                                         #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  #geom_point(aes(fill = sampletype, shape = filtertype), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23))+
  #scale_fill_manual(values = sample_colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))   

#import to metacoder

# Parse example dataset

mcdata <- parse_phyloseq(physeq)

# Plot data
heat_tree(x,
          node_size = n_obs,
          node_color = n_obs,
          node_label = taxon_names,
          tree_label = taxon_names)


## Fish only

phyfish = merge_phyloseq(subset_taxa(physeq, class=="Hyperoartia"), subset_taxa(physeq, class=="Actinopteri"))
phyfish = merge_samples(phyfish, "Site")

mcdata <- parse_phyloseq(phyfish)
heat_tree(mcdata,
          node_size = n_obs,
          node_color = n_obs,
          node_label = taxon_names,
          tree_label = taxon_names)

set.seed(1) # This makes the plot appear the same each time it is run 

heat_tree(mcdata, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = 'Big Stream')#, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

