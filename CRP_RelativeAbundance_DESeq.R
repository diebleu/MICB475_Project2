#load packages
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggforce)
library(DESeq2)
library(indicspecies)
library(BiocManager)
library(microbiome)
library(ggpubr)

anemia_metadata <- read_delim("anemia_metadata.txt", delim="\t")
tax <- read_delim("taxonomy.tsv", delim="\t")

#only keep healthy infant samples
filtered_meta_healthy <- filter(anemia_metadata, anemia == "normal", parasites == "N")
# filter out samples with no crp value
filtered_meta_healthy_crp <- filtered_meta_healthy[!is.na(as.numeric(filtered_meta_healthy$crp)),]

#make metadata file for only 6 month infant samples
meta_6m_crp <-  filter(filtered_meta_healthy_crp, age_months == 6)
#make metadata file for only month infant samples
meta_12m_crp <- filter(filtered_meta_healthy_crp, age_months == 12)

##create column based on clinical crp value - 6m
meta_6m_crp$crp_median <- ifelse(meta_6m_crp$crp >= 1,"Above", "Below")
#create column based on clinical crp value - 12m
meta_12m_crp$crp_median <- ifelse(meta_12m_crp$crp >= 1,"Above", "Below")

# Select columns
meta_6m_filt_crp <- select(meta_6m_crp, "#SampleID", "host_subject_id", "sex", "crp", "crp_status", "crp_median", "infection_status")
meta_12m_filt_crp <- select(meta_12m_crp, "#SampleID", "host_subject_id", "sex", "crp", "crp_status", "crp_median", "infection_status")

#save files
save(meta_6m_filt_crp, file = "crp_sorted_metadata_6m.RData")
save(meta_12m_filt_crp, file = "crp_sorted_metadata_12m.RData")

# Save as txt files - NOT WORKING
write.table(meta_12m_filt_crp, file = "crp_sorted_metadata_12m.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(meta_6m_filt_crp, file = "crp_sorted_metadata_6m.txt", sep="\t", quote=FALSE, row.names=FALSE)


#load in metadata and taxonomy tsvs
otu <- read_delim(file="table_250.tsv.txt", delim = "\t", skip=1)
tax <- read_delim(file = "taxonomy.tsv", delim="\t")
meta_6m <- read_delim("crp_sorted_metadata_6m.txt", delim="\t")
meta_12m <- read_delim("crp_sorted_metadata_12m.txt", delim="\t")
phylotree <- read.tree(file = "tree.nwk")

#rename first column of metadata
names(meta_6m)[names(meta_6m) == '#SampleID'] <- 'sampleid'
names(meta_12m)[names(meta_12m) == '#SampleID'] <- 'sampleid'

#filtering otu and tax files
otu_filt_6m <- otu %>% select("#OTU ID", one_of(meta_6m$sampleid))
otu_filt_12m <- otu %>% select("#OTU ID", one_of(meta_12m$sampleid))

tax_sep <- tax %>%
  separate(col=Taxon, sep=";", into = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species"))


#### Format OTU table ####

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df_6m <- as.data.frame(meta_6m[,-1])
samp_df_12m <- as.data.frame(meta_12m[,-1])

# Make sampleids the rownames
rownames(samp_df_6m)<- meta_6m$sampleid
rownames(samp_df_12m) <- meta_12m$sampleid

# Make phyloseq sample data with sample_data() function
SAMP_6m <- sample_data(samp_df_6m)
SAMP_12m <- sample_data(samp_df_12m)
class(SAMP_6m)
class(SAMP_12m)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence) %>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$"Feature ID"
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

class(phylotree)

#### Create phyloseq object ####
# Merge all into a phyloseq object
crp_infant_6m <- phyloseq(OTU, SAMP_6m, TAX, phylotree)
crp_infant_12m <- phyloseq(OTU, SAMP_12m, TAX, phylotree)




# Remove non-bacterial sequences, if any
crp_infant_6m_filt <- subset_taxa(crp_infant_6m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
crp_infant_12m_filt <- subset_taxa(crp_infant_12m,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Remove ASVs that have less than 5 counts total
crp_infant_6m_filt_nolow <- filter_taxa(crp_infant_6m_filt, function(x) sum(x)>5, prune = TRUE)
crp_infant_12m_filt_nolow <- filter_taxa(crp_infant_12m_filt, function(x) sum(x)>5, prune = TRUE)

# Genus Level
crp_infant_6m_genus <- tax_glom(crp_infant_6m_filt_nolow, "Genus")
crp_infant_12m_genus <- tax_glom(crp_infant_12m_filt_nolow, "Genus")

# Remove samples with less than 100 reads
crp_infant_6m_filt_nolow_samps <- prune_samples(sample_sums(crp_infant_6m_genus)>100, crp_infant_6m_filt_nolow)
crp_infant_12m_filt_nolow_samps <- prune_samples(sample_sums(crp_infant_12m_genus)>100, crp_infant_12m_filt_nolow)

# Remove samples where crp is na
crp_infant_6m_final <- subset_samples(crp_infant_6m_filt_nolow_samps, !is.na(crp) )
crp_infant_12m_final <- subset_samples(crp_infant_12m_filt_nolow_samps, !is.na(crp) )


################ Filtering Step 12months ###################

# Abundance filter: must be above 0.01% mean abundance
crp_infant_12m_final.rel = filter_taxa(crp_infant_12m_final, function(x) mean(x) > 1e-4, TRUE) 

# For each microbe, determine the proportion of samples it's present in
prevalencedf = apply(otu_table(crp_infant_12m_final.rel),MARGIN = 1,FUN = function(x)sum(x > 0)/nsamples(crp_infant_12m_final.rel))
prevalencedf = prevalencedf[prevalencedf>0.3] # Using 30% cutoff here

# Filter to only include these taxa
crp_infant_12m_final.rel = prune_taxa(names(prevalencedf), crp_infant_12m_final.rel)



##################### DESEQ ANALYSIS #############################

#setting random seed 
set.seed(1)

#Filtering Genus-level
# Genus Level
crp_infant_12m_genus <- tax_glom(crp_infant_12m_final.rel, "Genus")

# Remove samples with less than 100 reads
crp_infant_12m_filt_nolow_samps_DESeq <- prune_samples(sample_sums(crp_infant_12m_genus)>100, crp_infant_12m_genus)

# Remove samples where crp is na
crp_infant_12m_final_DESeq <- subset_samples(crp_infant_12m_filt_nolow_samps_DESeq, !is.na(crp) )

######### DESEQ ####################



#DESEq (No need to rarefaction [use infant_6m_final])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
infant_12m_plus1 <- transform_sample_counts(crp_infant_12m_genus, function(x) x+1)
### above could be infant_6m_plus1 <- transform_sample_counts(crp_infant_6m_final.rel, function(x) x+1)
infant_12m_deseq <- phyloseq_to_deseq2(infant_12m_plus1, ~crp_median)
DESEQ_infant_12m <- DESeq(infant_12m_deseq)

#Viewing Results
res_12m <- results(DESEQ_infant_12m, tidy=TRUE)
View(res_12m)

#Filtering out the ASV results with N.A. values for p adjusted values
filtered_res_12m <- drop_na(res_12m, padj)
view(filtered_res_12m)

#### Visualizing the results of DESEq ####
#Chose p value <0.05, and significant fold change >2

##Volcano Plot##
vol_plot_12m <- ggplot(filtered_res_12m, aes(x=log2FoldChange, y=-log(padj))) +
  geom_point()
vol_plot_sig_12m <- filtered_res_12m |>
  mutate(significant = padj<0.05 & abs(log2FoldChange) >2) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), col= significant)) +
  geom_point()

#Significant ASVs table
sigASVs_12m <- filtered_res_12m |>
  filter(padj <0.05 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

view(sigASVs_12m)

#BarPlot
#removed any Genus that has NA

sigASVs_vec_12m <- sigASVs_12m |>
  pull(ASV)

infant_12m_DESeq <- prune_taxa(sigASVs_vec_12m, crp_infant_12m_genus)
### or infant_6m_DESeq <- prune_taxa(sigASVs_vec_6m, crp_infant_6m_final_DESeq)

sigASVs_12m <- tax_table(infant_12m_DESeq) |>
  as.data.frame() |>
  rownames_to_column(var="ASV") |>
  right_join(sigASVs_12m) |>
  arrange(log2FoldChange) |>
  mutate(Genus = make.unique(Genus)) |>
  mutate(Genus = factor(Genus, levels= unique(Genus)))|>
  drop_na(Genus) |>
  filter(Genus != "NA.2" & Genus != "NA.1")

view(sigASVs_12m)

#crp_infant_12m_final.rel is the input for DESeq
#sigASVs_12m is the output from DESeq
#12 MONTHS
sigASVs_12m = sigASVs_12m %>% arrange(log2FoldChange) %>% mutate(Genus = factor(Genus,levels=.$Genus))

DESeq_crp_12m.compositional<- microbiome::transform(crp_infant_12m_final.rel, "compositional")

DESeq_crp_12m_psmelt <- psmelt(DESeq_crp_12m.compositional)

nrow(DESeq_crp_12m_psmelt)
ncol(DESeq_crp_12m_psmelt)
colnames(DESeq_crp_12m_psmelt)
view(DESeq_crp_12m_psmelt)

#renaming OTU column to ASV
colnames(DESeq_crp_12m_psmelt)[1] = "ASV"

DESeq_crp_12m_psmelt_sig <- DESeq_crp_12m_psmelt %>% filter(ASV %in% sigASVs_12m$ASV[sigASVs_12m$padj < 0.05])

nrow(DESeq_crp_12m_psmelt_sig)

#Adding Pseudocount
DESeq_crp_12m_psmelt_sig1 = DESeq_crp_12m_psmelt_sig %>% mutate(Abundance = Abundance + min(.$Abundance[.$Abundance>0]))

DESeq_crp_12m_psmelt_sig1 %>% 
  ggplot(aes(Genus,Abundance,fill=crp_median)) +
  geom_boxplot() +
  scale_y_log10(expand = expansion(mult = 0.1)) +
  theme_classic(base_size = 16) +
  scale_fill_discrete(labels=c('High', 'Low')) +
  labs(fill='CRP Level') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggpubr::stat_compare_means(aes(group=crp_median),label='p.signif',size=5)





################ Filtering Step 6months ###################

# Abundance filter: must be above 0.01% mean abundance
crp_infant_6m_final.rel = filter_taxa(crp_infant_6m_final, function(x) mean(x) > 1e-4, TRUE) 

# For each microbe, determine the proportion of samples it's present in
prevalencedf = apply(otu_table(crp_infant_6m_final.rel),MARGIN = 1,FUN = function(x)sum(x > 0)/nsamples(crp_infant_6m_final.rel))
prevalencedf = prevalencedf[prevalencedf>0.3] # Using 30% cutoff here

# Filter to only include these taxa
crp_infant_6m_final.rel = prune_taxa(names(prevalencedf), crp_infant_6m_final.rel)



##################### DESEQ ANALYSIS #############################

#setting random seed 
set.seed(1)

#Filtering Genus-level
# Genus Level
crp_infant_6m_genus <- tax_glom(crp_infant_6m_final.rel, "Genus")

# Remove samples with less than 100 reads
crp_infant_6m_filt_nolow_samps_DESeq <- prune_samples(sample_sums(crp_infant_6m_genus)>100, crp_infant_6m_genus)

# Remove samples where crp is na
crp_infant_6m_final_DESeq <- subset_samples(crp_infant_6m_filt_nolow_samps_DESeq, !is.na(crp) )

######### DESEQ ####################



#DESEq (No need to rarefaction [use infant_6m_final])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
infant_6m_plus1 <- transform_sample_counts(crp_infant_6m_genus, function(x) x+1)
### above could be infant_6m_plus1 <- transform_sample_counts(crp_infant_6m_final.rel, function(x) x+1)
infant_6m_deseq <- phyloseq_to_deseq2(infant_6m_plus1, ~crp_median)
DESEQ_infant_6m <- DESeq(infant_6m_deseq)

#Viewing Results
res_6m <- results(DESEQ_infant_6m, tidy=TRUE)

#Filtering out the ASV results with N.A. values for p adjusted values
filtered_res_6m <- drop_na(res_6m, padj)
view(filtered_res_6m)

#### Visualizing the results of DESEq ####
#Chose p value <0.05, and significant fold change >2

##Volcano Plot##
vol_plot_6m <- ggplot(filtered_res_6m, aes(x=log2FoldChange, y=-log(padj))) +
  geom_point()
vol_plot_sig_6m <- filtered_res_6m |>
  mutate(significant = padj<0.05 & abs(log2FoldChange) >2) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), col= significant)) +
  geom_point()

#Significant ASVs table
sigASVs_6m <- filtered_res_6m |>
  filter(padj <0.05 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

view(sigASVs_6m)

#BarPlot
#removed any Genus that has NA

sigASVs_vec_6m <- sigASVs_6m |>
  pull(ASV)

infant_6m_DESeq <- prune_taxa(sigASVs_vec_6m, crp_infant_6m_genus)
### or infant_6m_DESeq <- prune_taxa(sigASVs_vec_6m, crp_infant_6m_final_DESeq)

sigASVs_6m <- tax_table(infant_6m_DESeq) |>
  as.data.frame() |>
  rownames_to_column(var="ASV") |>
  right_join(sigASVs_6m) |>
  arrange(log2FoldChange) |>
  mutate(Genus = make.unique(Genus)) |>
  mutate(Genus = factor(Genus, levels= unique(Genus)))|>
  drop_na(Genus) |>
  filter(Genus != "NA.2" & Genus != "NA.1")

view(sigASVs_6m)

#crp_infant_6m_final.rel is the input for DESeq
#sigASVs_6m is the output from DESeq
#6 MONTHS
sigASVs_6m = sigASVs_6m %>% arrange(log2FoldChange) %>% mutate(Genus = factor(Genus,levels=.$Genus))

DESeq_crp_6m.compositional<- microbiome::transform(crp_infant_6m_final.rel, "compositional")

DESeq_crp_6m_psmelt <- psmelt(DESeq_crp_6m.compositional)

nrow(DESeq_crp_6m_psmelt)
ncol(DESeq_crp_6m_psmelt)
colnames(DESeq_crp_6m_psmelt)
view(DESeq_crp_6m_psmelt)

#renaming OTU column to ASV
colnames(DESeq_crp_6m_psmelt)[1] = "ASV"

DESeq_crp_6m_psmelt_sig <- DESeq_crp_6m_psmelt %>% filter(ASV %in% sigASVs_6m$ASV[sigASVs_6m$padj < 0.05])

nrow(DESeq_crp_6m_psmelt_sig)

#Adding Pseudocount
DESeq_crp_6m_psmelt_sig1 = DESeq_crp_6m_psmelt_sig %>% mutate(Abundance = Abundance + min(.$Abundance[.$Abundance>0]))

DESeq_crp_6m_psmelt_sig1 %>% 
  ggplot(aes(Genus,Abundance,fill=crp_median)) +
  geom_boxplot() +
  scale_y_log10(expand = expansion(mult = 0.1)) +
  theme_classic(base_size = 16) +
  scale_fill_discrete(labels=c('High', 'Low')) +
  labs(fill='CRP Level') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggpubr::stat_compare_means(aes(group=crp_median),label='p.signif',size=5)












