########### WORKSHOP DAY 4 #############
######## Further analysis & Figure Generation ############

#First install and load relevant libraries
install.packages("ggplot2")
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq", force = TRUE)

library(phyloseq)

BiocManager::install("decontam")

library(decontam)

library(devtools)

install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)

# We will now demonstrate how to use decontam to filter suspected contaminants
# You can do it one of two ways: by using your negative control OR library concentrations

#we will use decontam's pre-loaded dataset for this
ps <- readRDS(system.file("extdata", "MUClite.rds", package="decontam"))
ps

#Inspect the object; we have both the library quantity and controls
head(sample_data(ps))

#Inspect Library sizes first
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

#As you may expect, the true samples have many more reads than the control samples

####### Method 1: Frequency #########
#quant.reading is the column we want here
contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")
head(contamdf.freq)

# $p  containts the probability that was used for classifying contaminants
# TRUE = probable contaminant
#how many contaminants were identified by the analysis?

table(contamdf.freq$contaminant)

plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="quant_reading") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),3)], conc="quant_reading") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")


########## Method 2: Prevalence Method ##############
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

#This method found more contaminants than the frequency method
#We can change the threshold to see how it impacts contaminant identification

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

#Plot contams


#Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


#We will use the prevalence method to remove contaminants from the dataset
#Now need the same OTUs to disappear from the main physeq object
physeq_clean_decontam <- prune_taxa(contamdf.prev$contaminant == "FALSE", ps)
physeq_clean_decontam #check head(ps) to make sure OTU count decreased
ps

##################### Statistical Analyses with Phyloseq #######################

library(phyloseq)

physeq_clean_decontam #we will use the contaminant free object

#This is a big dataset, let's subset and only look at Tongue samples
tongue<-subset_samples(physeq_clean_decontam, Habitat=="Tongue")

#Convert to relative abundance
relative  = transform_sample_counts(tongue, function(OTU) OTU / sum(OTU))

#Visually inspect the top 50 OTUs
Top50OTUs = names(sort(taxa_sums(relative), TRUE)[1:50])

comparetop50 = prune_taxa(Top50OTUs, relative)

plot_bar(relative, fill = "Phylum", title = "Bacterial Phylum by Sample")
plot_bar(comparetop50, fill = "Phylum", title = "Bacterial Phylum by Sample")

#In this plot, each of the segments is 1 OTU. To make it prettier, we can collapse
#OTUs with the same taxonomy to a single segment in the bar plot.

plot_bar(comparetop50, fill = "Phylum", title = "Bacterial Phylum by Sample") +
  geom_bar(aes(), stat="identity", position="stack")

#We can even supply custom colors if we don't like the defaults
# viridis, wesanderson

values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darksalmon", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon")

plot_bar(comparetop50, fill = "Phylum", title = "Bacterial Phylum by Sample") +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = values)

#And separate plots into facets, in this case by sequencing plate
plot_bar(comparetop50, fill = "Phylum", title = "Bacterial Phylum by Sample") +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = values) +
  facet_grid(~PlateNumber)

##################### ALPHA DIVERSITY ####################################
#For this, do NOT use data that has already had singletons removed or some
#transformation to the data
require(wesanderson)

richness<-estimate_richness(physeq = physeq_clean_decontam, split = T, measures = c("Chao1","Shannon"))
richness$Plate_Number<- physeq_clean_decontam@sam_data$PlateNumber

Foxy<- wes_palette("FantasticFox1", 6, type = "continuous")
shannon.plot<-ggplot(richness, aes(Plate_Number, Shannon, fill= Plate_Number))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.5) +  #whiskers
  geom_boxplot(outlier.shape=1, outlier.colour = "red") + scale_fill_manual(values = Foxy) +
  theme_classic() + ggtitle("Plate Number")

#Which plate has the highest Shannon Richness?

###################### BETA DIVERSITY: ORDINATION & PERMANOVA #########################
#For this one, we want to use the transformed data (relative abundances)
#Set seed for reproducibility
set.seed(42)

rare<- rarefy_even_depth(physeq_clean_decontam)

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(rare))

# Calculate distance matrix
Dist.matrix <- distance(physeq_clean_decontam, method = "bray")

# PERMANOVA test- Habitat
require(vegan)
Habitat.Plate<-adonis2(formula = Dist.matrix ~ Habitat + PlateNumber, data = sampledf, permutations = 999)
Habitat
Habitat.Plate
#set strata = ColumnID to control for non-test variables

#Beta dispersion test: homogeneity of variance is an assumption of PERMANOVA
beta <- betadisper(Dist.matrix, sampledf$Habitat)
permutest(beta)

#Do we satisfy the assumption?

#Now we can plot the singnificant result
ord.bray <- ordinate(
  physeq = rare, 
  method = "PCoA", 
  distance = "bray"
)

betadivplot<-plot_ordination(
  physeq = rare,
  ordination = ord.bray,
  axes = c(1,2), 
  color = "Habitat")  +
  theme_classic() +
  geom_jitter()

betadivplot

##############Additional Utilities with MicrobiomeSeq ######################
require(microbiomeSeq)

#Load dataset on microbial communities in pit latrines
data(pitlatrine)

#We created an alpha diversity boxlpot in phyloseq, but maybe we'd like to statistically
# compare alpha diversity among groups using ANOVA

p<- plot_anova_diversity(pitlatrine, method = c("richness", "simpson", "shannon"), 
                         grouping_column = "Country", pValueCutoff = 0.05)
p


#Local contribution to beta diversity (LCBD)

physeq<- normalise_data(physeq = pitlatrine, norm.method = "relative") #also has more normalizing functions

p<- plot_taxa(physeq = physeq, grouping_column = "Country", method = "hellinger", 
              number.taxa = 21, filename = NULL)
p
