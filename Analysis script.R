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
library(ggrepel)

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
y <- taxonomy_table(x, use_ranks = c("superkingdom", "kingdom", "phylum", "class", "order",
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
#Fix names of waterways Kyeburn, Linnburn, Loganburn, Taiari, Waipouri
sampledata <- sampledata %>% 
  mutate(Waterway = str_replace(Waterway, "Kye burn", "Kyeburn")) %>% 
  mutate(Waterway = str_replace(Waterway, "Hog burn", "Hog Burn")) %>% 
  mutate(Waterway = str_replace(Waterway, "Linnburn river", "Linnburn")) %>%
  mutate(Waterway = str_replace(Waterway, "Logan Burn", "Loganburn")) %>%
  mutate(Waterway = str_replace(Waterway, "Taiari river", "Taiari River")) %>%
  mutate(Waterway = str_replace(Waterway, "Taiari River", "Taiari")) %>%
  mutate(Waterway = str_replace(Waterway, "Waipouri River", "Waipouri Stream")) %>%
  mutate(Waterway = str_replace(Waterway, "Waihola", "Waihora")) %>%
  mutate(Waterway = str_replace(Waterway, "Waipouri Stream", "Waipouri")) %>%
  mutate(Site = str_replace(Site, "Taiari River 1", "T1")) %>% 
  mutate(Site = str_replace(Site, "Taiari River 2", "T2")) %>%
  mutate(Site = str_replace(Site, "Taiari 3", "T3")) %>%
  mutate(Site = str_replace(Site, "Taiari 5", "T4")) %>%
  mutate(Site = str_replace(Site, "Taiari 6", "T5")) %>%
  mutate(Site = str_replace(Site, "Taiari 7", "T6")) %>%
  mutate(Site = str_replace(Site, "Taiari River 8", "T7")) %>%
  mutate(Site = str_replace(Site, "Taiari 9", "T8")) %>%
  mutate(Site = str_replace(Site, "Taiari 10", "T9")) %>%
  mutate(Site = str_replace(Site, "Taieri 10.5", "T10")) %>%
  mutate(Site = str_replace(Site, "Taiari River 11", "T11")) %>%
  mutate(Site = str_replace(Site, "Taiari 11.5", "T12")) %>%
  mutate(Site = str_replace(Site, "Taiari 12", "T13")) %>%
  mutate(Site = str_replace(Site, "Taiari 13", "T14")) %>%
  mutate(Site = str_replace(Site, "Taiari 14", "T15"))
  
#add to phyloseq object
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

all_pcoa <- ordinate(physeq, "PCoA", "bray")
p1 = plot_ordination(physeq, all_pcoa, type="taxa", color="phylum", title="Bacterial phyla")
print(p1)
p2 = plot_ordination(physeq, all_pcoa, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Waterway)) + geom_point(size=5) + ggtitle("samples")

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

### Bacteria

#Filter to only bacteria
bac <- subset_taxa(physeq, superkingdom == "Bacteria")
#Agglomerate bacteria taxa to family level
bac_glom <- tax_glom(bac, taxrank="family")
#remove 'kingdom' rank
tax_table(bac_glom) <- tax_table(bac_glom)[,c(1,3:6)]

#Plot barplots of relevant taxa ranks
plot_bar(bac_glom, x="Site", fill="phylum")

plot_bar(subset_taxa(bac_glom, phylum == "Proteobacteria"), x="Site", fill="order", facet_grid = ~class)
plot_bar(subset_taxa(bac_glom, phylum == "Chloroflexi"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Bacteroidetes"), x="Site", fill="order", facet_grid = ~class)
plot_bar(subset_taxa(bac_glom, phylum == "Spirochaetes"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Cyanobacteria"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Actinobacteria"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Deinococcus-Thermus"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Firmicutes"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Tenericutes"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Nitrospirae"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Acidobacteria"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Verrucomicrobia"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Chlamydiae"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Lentisphaerae"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Fusobacteria"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Planctomycetes"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Elusimicrobia"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Synergistetes"), x="Site", fill="class")
plot_bar(subset_taxa(bac_glom, phylum == "Armatimonadetes-Thermus"), x="Site", fill="class")

# Conduct ordination

bac.ord <- ordinate(bac_glom, "NMDS", "bray")
p1 = plot_ordination(bac_glom, bac.ord, type="taxa", color="phylum", title="Bacterial phyla")
print(p1)

p2 = plot_ordination(bac_glom, bac.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

proteo <- subset_taxa(bac_glom, phylum == "Proteobacteria")
proteo.ord <- ordinate(proteo, "NMDS", "bray")
p1 = plot_ordination(proteo, proteo.ord, type="taxa", color="order", title="Proteobacteria order")
print(p1)

p2 = plot_ordination(proteo, proteo.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

bacter <- subset_taxa(bac_glom, phylum == "Bacteroidetes")
bacter.ord <- ordinate(bacter, "NMDS", "bray")
p1 = plot_ordination(bacter, bacter.ord, type="taxa", color="order", title="Proteobacteria order")
print(p1)

p2 = plot_ordination(bacter, bacter.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

# Filter to only mainstem river samples
bac_tai <- subset_samples(bac, Waterway == "Taiari")
bac_tai <- phyloseq::subset_samples(bac_tai, Site != "T1")
bac_tai <- phyloseq::subset_samples(bac_tai, Site != "T2")
bac_tai <- phyloseq::subset_samples(bac_tai, Site != "T15")# Remove Taiari River 1 as it distorts the data
bac.tai.ord <- ordinate(bac_tai , "PCoA", "bray")
p2 <- plot_ordination(bac_tai, bac.tai.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

bac.taiari <- data.frame(cbind(bac.tai.ord$vectors[,1:2], data.frame(sample_data(bac_tai))$Site))
colnames(bac.taiari) <- c("Ax1", "Ax2", "Site")
bac.taiari$Ax1 <- as.numeric(bac.taiari$Ax1)
bac.taiari$Ax2 <- as.numeric(bac.taiari$Ax2)
bac.tai.plot.data <- cbind(aggregate(x = bac.taiari$Ax1, by = list(bac.taiari$Site), FUN = "mean"),
                           aggregate(x = bac.taiari$Ax2, by = list(bac.taiari$Site), FUN = "mean"))
bac.tai.plot.data <- as.data.frame(cbind(bac.tai.plot.data[,2], bac.tai.plot.data[,4], bac.tai.plot.data[,1]))
colnames(bac.tai.plot.data) <- c("Ax1", "Ax2", "Site")
bac.tai.plot.data$Ax1 <- as.numeric(bac.tai.plot.data$Ax1)
bac.tai.plot.data$Ax2 <- as.numeric(bac.tai.plot.data$Ax2)
bac.tai.plot.data$Site <- as.factor(bac.tai.plot.data$Site)
bac.tai.plot.data <- arrange(bac.tai.plot.data, factor(Site, levels = c("T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13", "T14")))

bac.t.ord <- bac.tai.plot.data %>%
  ggplot() +
  theme_bw() +
  geom_point(aes(x=Ax1, y=Ax2, colour=Site), size =3, show.legend = TRUE) +
  geom_path(aes(x=Ax1, y=Ax2), linetype="dotted") +
  geom_text(aes(x=Ax1, y=Ax2),
            label = bac.tai.plot.data$Site,
            nudge_x = 0.01, nudge_y = 0.01,
            check_overlap = F)


### Fish


#Filter to only fish
fish <- merge_phyloseq(subset_taxa(physeq, class=="Hyperoartia"), subset_taxa(physeq, class=="Actinopteri"))

#remove 'superkingdom' nad 'kingdom' ranks
tax_table(fish) <- tax_table(fish)[,c(3:8)]

#Plot barplots of relevant taxa ranks
plot_bar(fish, x="Site", fill="family") #, facet_grid = ~order)

plot_bar(subset_taxa(fish, family == "Anguillidae"), x="Site", fill="species")
plot_bar(subset_taxa(fish, family == "Geotriidae"), x="Site", fill="species")
plot_bar(subset_taxa(fish, family == "Percidae"), x="Site", fill="species")
plot_bar(subset_taxa(fish, family == "Scombridae"), x="Site", fill="genus") #Contamination 3 o'clock
plot_bar(subset_taxa(fish, family == "Latridae"), x="Site", fill="species") # Blue Moki T1
plot_bar(subset_taxa(fish, family == "Retropinnidae"), x="Site", fill="species") # Smelt T2
plot_bar(subset_taxa(fish, family == "Tripterygiidae"), x="Site", fill="species") # Triplefin T1
plot_bar(subset_taxa(fish, family == "Cheimarrichthyidae"), x="Site", fill="species") # Torrentfish 3 oclock, Meggat Burn and T3


plot_bar(subset_taxa(fish, family == "Cheimarrichthyidae"), x="Site", fill="species") # Torrentfish 3 oclock, Meggat Burn and T3

plot_bar(
  merge_samples(prune_taxa(
  subset_samples(fish, Site == "T3"), "Site"), taxa_sums > 0),  fill="species", 
  facet_grid = ~species)


# Conduct ordination
fish.0 <- prune_samples(sample_sums(fish) > 0, fish)
fish.ord <- ordinate(fish.0, "NMDS", "bray")
p1 = plot_ordination(fish.0, fish.ord, type="taxa", color="species", title="Fish species")
print(p1)
p2 = plot_ordination(fish.0, fish.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

galax <- subset_taxa(fish, genus == "Galaxias")
galax.0 <- prune_samples(sample_sums(galax) > 0, galax)
galax.ord <- ordinate(galax.0, "NMDS", "bray")
p1 = plot_ordination(galax.0, galax.ord, type="taxa", color="species", title="Galaxias species")
print(p1)
p2 = plot_ordination(galax.0, galax.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

gobio <- subset_taxa(fish, genus == "Gobiomorphus")
gobio.0 <- prune_samples(sample_sums(gobio) > 0, gobio)
gobio.ord <- ordinate(gobio.0, "NMDS", "bray")
p1 = plot_ordination(gobio.0, gobio.ord, type="taxa", color="species", title="Gobiomorphus species")
print(p1)
p2 = plot_ordination(gobio.0, gobio.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

heat_tree(parse_phyloseq(fish),
          node_size = n_obs,
          node_color = n_obs,
          node_label = taxon_names,
          tree_label = taxon_names)

### Diatoms

#Filter to only diatoms
diatom <- subset_taxa(physeq, phylum=="Bacillariophyta")

#remove 'superkingdom' and 'kingdom' ranks
tax_table(diatom) <- tax_table(diatom)[,c(3:8)]

#Plot barplots of relevant taxa ranks
plot_bar(diatom, x="Site", fill="order")#, facet_grid = ~order)

# Conduct ordination
diatom.0 <- prune_samples(sample_sums(diatoms) > 0, diatoms)
diatom.ord <- ordinate(diatom.0, "NMDS", "bray")
p1 = plot_ordination(diatom.0, diatom.ord, type="taxa", color="family", title="Diatom family")
print(p1)
p2 = plot_ordination(diatom.0, diatom.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

#Filter to only ciliates
ciliate <- subset_taxa(physeq, phylum=="Ciliophora")
ciliate <- tax_glom(ciliate, taxrank="genus")

#remove 'superkingdom' and 'kingdom' ranks
tax_table(ciliate) <- tax_table(ciliate)[,c(3:7)]

#Plot barplots of relevant taxa ranks
plot_bar(ciliate, x="Site", fill="order")#, facet_grid = ~order)

# Conduct ordination
ciliate.0 <- prune_samples(sample_sums(ciliate) > 0, ciliate)
ciliate.ord <- ordinate(ciliate.0, "NMDS", "bray")
p1 = plot_ordination(ciliate.0, ciliate.ord, type="taxa", color="family", title="Giliate families")
print(p1)
p2 = plot_ordination(ciliate.0, ciliate.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

### Chlorophytes - Green Algae

#Filter to only green algae
chlorophytes <- subset_taxa(physeq, phylum=="Chlorophyta")
#chlorophytes <- tax_glom(chlorophytes, taxrank="genus")

#remove 'superkingdom' and 'kingdom' ranks
tax_table(chlorophytes) <- tax_table(chlorophytes)[,c(3:8)]

#Plot barplots of relevant taxa ranks
plot_bar(chlorophytes, x="Site", fill="order")#, facet_grid = ~order)

# Conduct ordination
chlorophytes.0 <- prune_samples(sample_sums(chlorophytes) > 0, chlorophytes)
chlorophytes.ord <- ordinate(chlorophytes.0, "NMDS", "bray")
p1 = plot_ordination(chlorophytes.0, chlorophytes.ord, type="taxa", color="family", title="Giliate families")
print(p1)
p2 = plot_ordination(chlorophytes.0, chlorophytes.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

#Filter to only insects
insects <- subset_taxa(physeq, class=="Insecta")
#chlorophytes <- tax_glom(chlorophytes, taxrank="genus")

#remove 'superkingdom' and 'kingdom' ranks
tax_table(insects) <- tax_table(insects)[,c(4:8)]

#Plot barplots of relevant taxa ranks
plot_bar(insects, x="Site", fill="order")#, facet_grid = ~order)

# Conduct ordination
insects.0 <- prune_samples(sample_sums(insects) > 0, insects)
insects.ord <- ordinate(insects.0, "NMDS", "bray")
p1 = plot_ordination(insects.0, insects.ord, type="taxa", color="order", title="Insects order")
print(p1)
p2 = plot_ordination(insects.0, insects.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

mcdata <- parse_phyloseq(insects)
heat_tree(mcdata,
          node_size = n_obs,
          node_color = n_obs,
          node_label = taxon_names,
          tree_label = taxon_names)


