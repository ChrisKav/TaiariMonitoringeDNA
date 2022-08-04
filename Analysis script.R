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

otumat <- data.frame(taxa_data[c(3,6:ncol(taxa_data))])

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

### SETUP TO IMPORT INTO PHYSLOSEQ
row.names(otumat) <- otumat$TaxID
otumat$TaxID <- NULL
otumat <- as.matrix(otumat)
OTU = otu_table(otumat, taxa_are_rows = TRUE)

### USE library(taxa and metacoder) to do taxonomy for phyloseq
x <- lookup_tax_data(taxdf$ID, type = "taxon_id") #takes a while to run
y <- taxonomy_table(x, use_ranks = c("superkingdom", "kingdom", "phylum", "class", "order",
        "family", "genus", "species"), add_id_col = TRUE)
y <- data.frame(y)
rownames(y) <- y$taxon_id
y$taxon_id <- NULL
taxmat <- as.matrix(y)

# write.csv(taxmat, "tax_table.csv") # Check

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

#saveRDS(physeq, file = "April 2022 data.rds")

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
plot_bar(subset_taxa(fish, family == "Salmonidae"), x="Site", fill="genus") # Brown trout everywhere!
plot_bar(subset_taxa(fish, family == "Galaxiidae"), x="Site", fill="species") # Galaxiids
plot_bar(subset_taxa(fish, family == "Eleotridae"), x="Site", fill="species") # Bullies
plot_bar(subset_taxa(fish, family == "Mugilidae"), x="Site", fill="species") # Yellow eyed mullet - T1 and T2, Waihola and Contour Channel
plot_bar(subset_taxa(fish, family == "Moridae"), x="Site", fill="genus") # T1 Rock cod
plot_bar(subset_taxa(fish, family == "Gempylidae"), x="Site", fill="species") # T1 Snoek
plot_bar(subset_taxa(fish, family == "Monacanthidae"), x="Site", fill="species") #T1 Leatherjacket
plot_bar(subset_taxa(fish, family == "Labridae"), x="Site", fill="species") #parrot fish T1
plot_bar(subset_taxa(fish, family == "Rhombosoleidae"), x="Site", fill="species") # flounders T1 and T3
plot_bar(subset_taxa(fish, family == "Hemiramphidae"), x="Site", fill="species") #T1 garfish
plot_bar(subset_taxa(fish, family == "Carangidae"), x="Site", fill="family") #Jacks

plot_bar(subset_taxa(fish, species == "Galaxias anomalus"), x="Site", fill="species") #CO Rundheads
plot_bar(subset_taxa(fish, species == "Galaxias argenteus"), x="Site", fill="species") #GKs
plot_bar(subset_taxa(fish, species == "Galaxias brevipinnis"), x="Site", fill="species") #Koaro
plot_bar(subset_taxa(fish, species == "Galaxias depressiceps"), x="Site", fill="species") # T Flatheads
plot_bar(subset_taxa(fish, species == "Galaxias eldoni"), x="Site", fill="species") #Eldons
plot_bar(subset_taxa(fish, species == "Galaxias fasciatus"), x="Site", fill="species") #Banded Kokopu
plot_bar(subset_taxa(fish, species == "Galaxias maculatus"), x="Site", fill="species") #Inaka
plot_bar(subset_taxa(fish, species == "Galaxias pullus"), x="Site", fill="species") #Duskys
plot_bar(subset_taxa(fish, species == "Galaxias sp. 'teviot'"), x="Site", fill="species") #Teviots - incorrect - anomalus
plot_bar(subset_taxa(fish, species == "Galaxias sp. D (Allibone et al., 1996)"), x="Site", fill="species") #Clutha flats

T1.fish <- subset_samples(fish, Site == c("T1"))
T1.fish <- prune_taxa(taxa_sums(T1.fish) >0, T1.fish)
T1.fish <- merge_samples(T1.fish, "Site")
plot_bar(T1.fish, fill="species", facet_grid = ~species)

T2.fish <- subset_samples(fish, Site == "T2")
T2.fish <- prune_taxa(taxa_sums(T2.fish) >0, T2.fish)
T2.fish <- merge_samples(T2.fish, "Site")
plot_bar(T2.fish, fill="species", facet_grid = ~species)

T3.fish <- subset_samples(fish, Site == "T3")
T3.fish <- prune_taxa(taxa_sums(T3.fish) >0, T3.fish)
T3.fish <- merge_samples(T3.fish, "Site")
plot_bar(T3.fish, fill="species", facet_grid = ~species)

T4.fish <- subset_samples(fish, Site == "T4")
T4.fish <- prune_taxa(taxa_sums(T4.fish) >0, T4.fish)
T4.fish <- merge_samples(T4.fish, "Site")
plot_bar(T4.fish, fill="species", facet_grid = ~species)

T5.fish <- subset_samples(fish, Site == "T5")
T5.fish <- prune_taxa(taxa_sums(T5.fish) >0, T5.fish)
T5.fish <- merge_samples(T5.fish, "Site")
plot_bar(T5.fish, fill="species", facet_grid = ~species)

T6.fish <- subset_samples(fish, Site == "T6")
T6.fish <- prune_taxa(taxa_sums(T6.fish) >0, T6.fish)
T6.fish <- merge_samples(T6.fish, "Site")
plot_bar(T6.fish, fill="species", facet_grid = ~species)

T7.fish <- subset_samples(fish, Site == "T7")
T7.fish <- prune_taxa(taxa_sums(T7.fish) >0, T7.fish)
T7.fish <- merge_samples(T7.fish, "Site")
plot_bar(T7.fish, fill="species", facet_grid = ~species)

T8.fish <- subset_samples(fish, Site == "T8")
T8.fish <- prune_taxa(taxa_sums(T8.fish) >0, T8.fish)
T8.fish <- merge_samples(T8.fish, "Site")
plot_bar(T8.fish, fill="species", facet_grid = ~species)

T9.fish <- subset_samples(fish, Site == "T9")
T9.fish <- prune_taxa(taxa_sums(T9.fish) >0, T9.fish)
T9.fish <- merge_samples(T9.fish, "Site")
plot_bar(T9.fish, fill="species", facet_grid = ~species)

T10.fish <- subset_samples(fish, Site == "T10")
T10.fish <- prune_taxa(taxa_sums(T10.fish) >0, T10.fish)
T10.fish <- merge_samples(T10.fish, "Site")
plot_bar(T10.fish, fill="species", facet_grid = ~species)

T11.fish <- subset_samples(fish, Site == "T11")
T11.fish <- prune_taxa(taxa_sums(T11.fish) >0, T11.fish)
T11.fish <- merge_samples(T11.fish, "Site")
plot_bar(T11.fish, fill="species", facet_grid = ~species)

T12.fish <- subset_samples(fish, Site == "T12")
T12.fish <- prune_taxa(taxa_sums(T12.fish) >0, T12.fish)
T12.fish <- merge_samples(T12.fish, "Site")
plot_bar(T12.fish, fill="species", facet_grid = ~species)

T13.fish <- subset_samples(fish, Site == "T13")
T13.fish <- prune_taxa(taxa_sums(T13.fish) >0, T13.fish)
T13.fish <- merge_samples(T13.fish, "Site")
plot_bar(T13.fish, fill="species", facet_grid = ~species)

T14.fish <- subset_samples(fish, Site == "T14")
T14.fish <- prune_taxa(taxa_sums(T14.fish) >0, T14.fish)
T14.fish <- merge_samples(T14.fish, "Site")
plot_bar(T14.fish, fill="species", facet_grid = ~species)

T15.fish <- subset_samples(fish, Site == "T15")
T15.fish <- prune_taxa(taxa_sums(T15.fish) >0, T15.fish)
T15.fish <- merge_samples(T15.fish, "Site")
plot_bar(T15.fish, fill="species", facet_grid = ~species)

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
          tree_label = taxon_names,
          layout = "davidson-harel")

### Molluscs
plot_bar(subset_taxa(physeq, family == "Lymnaeidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Planorbidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, species == "Physella acuta"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, genus == "Potamopyrgus"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, family == "Sphaeriidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, species == "Echyridella menziesii"), x="Site", fill="species")

### Diatoms

#Filter to only diatoms
diatom <- subset_taxa(physeq, phylum=="Bacillariophyta")

#remove 'superkingdom' and 'kingdom' ranks
tax_table(diatom) <- tax_table(diatom)[,c(3:8)]

#Plot barplots of relevant taxa ranks
plot_bar(diatom, x="Site", fill="order")#, facet_grid = ~order)

# Conduct ordination
diatom.0 <- prune_samples(sample_sums(diatom) > 0, diatom)
diatom.ord <- ordinate(diatom.0, "PCoA", "bray")
p1 = plot_ordination(diatom.0, diatom.ord, type="taxa", color="family", title="Diatom family")
print(p1)
p2 = plot_ordination(diatom.0, diatom.ord, type="samples", color="Site") 
p2 + geom_polygon(aes(fill=Site)) + geom_point(size=5) + ggtitle("samples")

plot_bar(subset_taxa(diatom, phylum == "Bacillariophyta"), x="Site", fill="family")
plot_bar(subset_taxa(diatom, class == "Bacillariophyceae"), x="Site", fill="family")
plot_bar(subset_taxa(diatom, family == "Gomphonemataceae"), x="Site", fill="genus")
plot_bar(subset_taxa(diatom, genus == "Didymosphenia"), x="Site", fill="species")
plot_bar(subset_taxa(diatom, class == "Coscinodiscophyceae"), x="Site", fill="order")
plot_bar(subset_taxa(diatom, family == "Cymbellaceae"), x="Site", fill="genus")
plot_bar(subset_taxa(diatom, phylum == "Bacillariophyta"), x="Site", fill="class")
plot_bar(subset_taxa(diatom, order == "Melosirales"), x="Site", fill="family")


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
          tree_label = taxon_names,
          layout = "davidson-harel")

### Arthropods
plot_bar(subset_taxa(physeq, class == "Branchiopoda"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Daphniidae"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, family == "Bosminidae"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, family == "Chydoridae"), x="Site", fill="species")

plot_bar(subset_taxa(physeq, class == "Ostracoda"), x="Site", fill="genus")

plot_bar(subset_taxa(physeq, class == "Malacostraca"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, species == "Paranephrops zealandicus"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, order == "Amphipoda"), x="Site", fill="genus")

plot_bar(subset_taxa(physeq, class == "Collembola"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, class == "Collembola"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, class == "Arachnida"), x="Site", fill="family")
plot_bar(subset_taxa(physeq, class == "Pycnogonida"), x="Site", fill="genus")

plot_bar(subset_taxa(physeq, class == "Insecta"), x="Site", fill="order")
plot_bar(subset_taxa(physeq, order == "Plecoptera"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, genus == "Sigara"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Coleoptera"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Trichoptera"), x="Site", fill="family")
plot_bar(subset_taxa(physeq, family == "Hydrobiosidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, genus == "Psilochorema"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, genus == "Costachorema"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, genus == "Hydrobiosis"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, family == "Leptoceridae"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, family == "Conoesucidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Hydroptilidae"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, family == "Hydropsychidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Philopotamidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Polycentropodidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Helicopsychidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Helicophidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Oeconesidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Diptera"), x="Site", fill="family")
plot_bar(subset_taxa(physeq, family == "Chironomidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Syrphidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Simuliidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Culicidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Ceratopogonidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Cecidomyiidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Psychodidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Limoniidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Ephydridae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Lonchopteridae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Syrphidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Megaloptera"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Ephemeroptera"), x="Site", fill="family")
plot_bar(subset_taxa(physeq, family == "Leptophlebiidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, genus == "Deleatidium"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, family == "Coloburiscidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Nesameletidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Ameletopsidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Oniscigastridae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, genus == "Neozephlebia"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, genus == "Zephlebia"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, genus == "Atalophlebioides"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Odonata"), x="Site", fill="genus")

### Cnidarians
plot_bar(subset_taxa(physeq, phylum == "Cnidaria"), x="Site", fill="order")
plot_bar(subset_taxa(physeq, class == "Myxozoa"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, class == "Hydrozoa"), x="Site", fill="order")
plot_bar(subset_taxa(physeq, order == "Anthoathecata"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, genus == "Hydra"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, order == "Limnomedusae"), x="Site", fill="genus")


### Proifera
plot_bar(subset_taxa(physeq, phylum == "Porifera"), x="Site", fill="order")

### Bryzoa
plot_bar(subset_taxa(physeq, phylum == "Bryozoa"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, genus == "Electra"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, family == "Plumatellidae"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, family == "Fredericellidae"), x="Site", fill="genus")

### Rotifera
plot_bar(subset_taxa(physeq, phylum == "Rotifera"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Ploima"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Flosculariaceae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Philodinida"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Adinetida"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Philodinavidae"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, genus == "Rotaria"), x="Site", fill="species")

### Platyhelminthes
plot_bar(subset_taxa(physeq, phylum == "Platyhelminthes"), x="Site", fill="family")
plot_bar(subset_taxa(physeq, class == "Rhabditophora"), x="Site", fill="family")
plot_bar(subset_taxa(physeq, order == "Tricladida"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Rhabdocoela"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Macrostomida"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Lecithoepitheliata"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, class == "Trematoda"), x="Site", fill="family")
plot_bar(subset_taxa(physeq, class == "Catenulida"), x="Site", fill="family")

plot_bar(subset_taxa(physeq, phylum == "Nemertea"), x="Site", fill="species")
plot_bar(subset_taxa(physeq, phylum == "Tardigrada"), x="Site", fill="family")

plot_bar(subset_taxa(physeq, phylum == "Annelida"), x="Site", fill="family")
plot_bar(subset_taxa(physeq, class == "Polychaeta"), x="Site", fill="order")
plot_bar(subset_taxa(physeq, order == "Tubificida"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Enchytraeida"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Lumbriculida"), x="Site", fill="genus")
plot_bar(subset_taxa(physeq, order == "Rhynchobdellida"), x="Site", fill="genus")

plot_bar(subset_taxa(physeq, phylum == "Gastrotricha"), x="Site", fill="species")

plot_bar(subset_taxa(physeq, phylum == "Nematoda"), x="Site", fill="family")

plot_bar(subset_taxa(physeq, class == "Amphibia"), x="Site", fill="genus") # Frogs at Eden Creek

###Birds

plot_bar(subset_taxa(physeq, genus == "Charadrius"), x="Site", fill="species") # New Zealand dotterel Kyeburn
plot_bar(subset_taxa(physeq, genus == "Hemiphaga"), x="Site", fill="species") # Kereru Silverstream and Traquair
plot_bar(subset_taxa(physeq, genus == "Zosterops"), x="Site", fill="species") # Silvereye
plot_bar(subset_taxa(physeq, genus == "Gerygone"), x="Site", fill="species") # Riroriro
plot_bar(subset_taxa(physeq, genus == "Anthornis"), x="Site", fill="species") # Korimako
plot_bar(subset_taxa(physeq, genus == "Petroica"), x="Site", fill="species") # Tomtit
plot_bar(subset_taxa(physeq, genus == "Phalacrocorax"), x="Site", fill="species") # Cormorant
plot_bar(subset_taxa(physeq, genus == "Microcarbo"), x="Site", fill="species") # Cormorant
plot_bar(subset_taxa(physeq, genus == "Egretta"), x="Site", fill="species") # White Heron
plot_bar(subset_taxa(physeq, genus == "Eudyptula"), x="Site", fill="species") # Little penguin at Waipouri
plot_bar(subset_taxa(physeq, genus == "Porphyrio"), x="Site", fill="species") # Pukeko
plot_bar(subset_taxa(physeq, genus == "Circus"), x="Site", fill="species") # Swamp Harrier

################################################################
library(microeco)
library(file2meco)
library(pheatmap)
library(magrittr)

#First remove kingdom rank and rename superkingdom to kingdom
bac2 <-  merge_phyloseq(subset_taxa(physeq, superkingdom=="Bacteria"), subset_taxa(physeq, superkingdom=="Archaea"))

tax_table(bac2) <- tax_table(bac2)[,c(1,3:8)]
colnames(tax_table(bac2)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")

heat_tree(parse_phyloseq(bac2),
          node_size = n_obs,
          node_color = n_obs,
          node_label = taxon_names,
          tree_label = taxon_names,
          layout = "davidson-harel")

#Transfer to meco
meco_dataset <- phyloseq2meco(bac2)
# create object of trans_func
t2 <- trans_func$new(meco_dataset)
# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database
t2$cal_spe_func(prok_database = "FAPROTAX")
# return t2$res_spe_func, 1 represent trait exists, 0 represent no or cannot confirmed.
t2$res_spe_func[1:5, 1:2]
# calculate the percentages for communities
# here do not consider the abundance
t2$cal_spe_func_perc(abundance_weighted = FALSE)
## The result table is stored in object$res_spe_func_perc ...
t2$res_spe_func_perc[1:5, 1:2]
# construct a network for the example
network <- trans_network$new(dataset = t2, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
network$cal_network(p_thres = 0.01, COR_cut = 0.7)
network$cal_module()
# convert module info to microtable object
meco_module <- network$trans_comm(use_col = "module")
meco_module_func <- trans_func$new(meco_module)
meco_module_func$cal_spe_func(prok_database = "FAPROTAX")
meco_module_func$cal_spe_func_perc(abundance_weighted = FALSE)
meco_module_func$plot_spe_func_perc()
# M represents module, ordered by the nodes number from high to low
# If you want to change the group list, reset the list t2$func_group_list
t2$func_group_list
# use show_prok_func to see the detailed information of prokaryotic traits
t2$show_prok_func("nitrate_reduction")
# then we try to correlate the res_spe_func_perc of communities to environmental variables
t3 <- trans_env$new(dataset = meco_dataset, add_data = data.frame(sample_data(bac2)[,c(1:3,6)]))
t3$cal_cor(add_abund_table = t2$res_spe_func_perc, cor_method = "spearman")
t3$plot_cor(pheatmap = TRUE)

View(t2$res_spe_func_perc)

m1 <- trans_abund$new(dataset = meco_dataset, taxrank = 'Phylum', groupmean = "Site", ntaxa = 5)
m1$plot_pie(facet_nrow = 5)





m1 <- trans_func$new(meco_dataset)
m1
m1$cal_tax4fun(folderReferenceData = "./SILVA123")

# must transpose to taxa row, sample column
pathway_file <- m1$tax4fun_path$Tax4FunProfile %>% t %>% as.data.frame
# filter rownames, only keep ko+number
rownames(pathway_file) %<>% gsub("(^.*);\\s.*", "\\1", .)
# load the pathway hierarchical metadata
data(Tax4Fun2_KEGG)
# further create a microtable object, familiar?
func1 <- microtable$new(otu_table = pathway_file, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = m1$sample_table)
print(func1)

func1$tidy_dataset()
# calculate abundance automatically at three levels: Level.1, Level.2, Level.3
func1$cal_abund()
print(func1)

# bar plot at Level.1
func2 <- trans_abund$new(func1, taxrank = "Level.3", groupmean = "Site")
func2$plot_bar(legend_text_italic = FALSE)

func2 <- trans_diff$new(dataset = func1, method = "lefse", group = "Site", alpha = 0.05, lefse_subgroup = NULL)
func2$plot_diff_bar(threshold = 4, width = 0.8)
