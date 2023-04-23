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


################ Filtering Step 12months ###################

# Abundance filter: must be above 0.01% mean abundance
agp_infant_12m_final.rel = filter_taxa(agp_infant_12m_final, function(x) mean(x) > 1e-4, TRUE) 

# For each microbe, determine the proportion of samples it's present in
prevalencedf = apply(otu_table(agp_infant_12m_final.rel),MARGIN = 1,FUN = function(x)sum(x > 0)/nsamples(agp_infant_12m_final.rel))
prevalencedf = prevalencedf[prevalencedf>0.3] # Using 30% cutoff here

# Filter to only include these taxa
agp_infant_12m_final.rel = prune_taxa(names(prevalencedf), agp_infant_12m_final.rel)



##################### DESEQ ANALYSIS #############################

#setting random seed 
set.seed(1)

#Filtering Genus-level
# Genus Level
agp_infant_12m_genus <- tax_glom(agp_infant_12m_final.rel, "Genus")

# Remove samples with less than 100 reads
agp_infant_12m_filt_nolow_samps_DESeq <- prune_samples(sample_sums(agp_infant_12m_genus)>100, agp_infant_12m_genus)

# Remove samples where agp is na
agp_infant_12m_final_DESeq <- subset_samples(agp_infant_12m_filt_nolow_samps_DESeq, !is.na(agp) )

######### DESEQ ####################



#DESEq (No need to rarefaction [use infant_6m_final])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
infant_12m_plus1 <- transform_sample_counts(agp_infant_12m_genus, function(x) x+1)
### above could be infant_6m_plus1 <- transform_sample_counts(agp_infant_6m_final.rel, function(x) x+1)
infant_12m_deseq <- phyloseq_to_deseq2(infant_12m_plus1, ~agp_clin)
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

infant_12m_DESeq <- prune_taxa(sigASVs_vec_12m, agp_infant_12m_genus)
### or infant_6m_DESeq <- prune_taxa(sigASVs_vec_6m, agp_infant_6m_final_DESeq)

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

#agp_infant_12m_final.rel is the input for DESeq
#sigASVs_12m is the output from DESeq
#12 MONTHS
sigASVs_12m = sigASVs_12m %>% arrange(log2FoldChange) %>% mutate(Genus = factor(Genus,levels=.$Genus))

DESeq_agp_12m.compositional<- microbiome::transform(agp_infant_12m_final.rel, "compositional")

DESeq_agp_12m_psmelt <- psmelt(DESeq_agp_12m.compositional)

nrow(DESeq_agp_12m_psmelt)
ncol(DESeq_agp_12m_psmelt)
colnames(DESeq_agp_12m_psmelt)
view(DESeq_agp_12m_psmelt)

#renaming OTU column to ASV
colnames(DESeq_agp_12m_psmelt)[1] = "ASV"

DESeq_agp_12m_psmelt_sig <- DESeq_agp_12m_psmelt %>% filter(ASV %in% sigASVs_12m$ASV[sigASVs_12m$padj < 0.05])

nrow(DESeq_agp_12m_psmelt_sig)

#Adding Pseudocount
DESeq_agp_12m_psmelt_sig1 = DESeq_agp_12m_psmelt_sig %>% mutate(Abundance = Abundance + min(.$Abundance[.$Abundance>0]))

DESeq_agp_12m_psmelt_sig1 %>% 
  ggplot(aes(Genus,Abundance,fill=agp_clin)) +
  geom_boxplot() +
  scale_y_log10(expand = expansion(mult = 0.1)) +
  theme_classic(base_size = 16) +
  scale_fill_discrete(labels=c('High', 'Low')) +
  labs(fill='AGP Level') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggpubr::stat_compare_means(aes(group=agp_clin),label='p.signif',size=5)
 




################ Filtering Step 6months ###################

# Abundance filter: must be above 0.01% mean abundance
agp_infant_6m_final.rel = filter_taxa(agp_infant_6m_final, function(x) mean(x) > 1e-4, TRUE) 

# For each microbe, determine the proportion of samples it's present in
prevalencedf = apply(otu_table(agp_infant_6m_final.rel),MARGIN = 1,FUN = function(x)sum(x > 0)/nsamples(agp_infant_6m_final.rel))
prevalencedf = prevalencedf[prevalencedf>0.3] # Using 30% cutoff here

# Filter to only include these taxa
agp_infant_6m_final.rel = prune_taxa(names(prevalencedf), agp_infant_6m_final.rel)



##################### DESEQ ANALYSIS #############################

#setting random seed 
set.seed(1)

#Filtering Genus-level
# Genus Level
agp_infant_6m_genus <- tax_glom(agp_infant_6m_final.rel, "Genus")

# Remove samples with less than 100 reads
agp_infant_6m_filt_nolow_samps_DESeq <- prune_samples(sample_sums(agp_infant_6m_genus)>100, agp_infant_6m_genus)

# Remove samples where agp is na
agp_infant_6m_final_DESeq <- subset_samples(agp_infant_6m_filt_nolow_samps_DESeq, !is.na(agp) )

######### DESEQ ####################



#DESEq (No need to rarefaction [use infant_6m_final])
#ERROR Message regarding genes having 0 reads.
#Need to add '1' read count to all genes [LIMITATION]
infant_6m_plus1 <- transform_sample_counts(agp_infant_6m_genus, function(x) x+1)
### above could be infant_6m_plus1 <- transform_sample_counts(agp_infant_6m_final.rel, function(x) x+1)
infant_6m_deseq <- phyloseq_to_deseq2(infant_6m_plus1, ~agp_clin)
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

infant_6m_DESeq <- prune_taxa(sigASVs_vec_6m, agp_infant_6m_genus)
### or infant_6m_DESeq <- prune_taxa(sigASVs_vec_6m, agp_infant_6m_final_DESeq)

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

#agp_infant_6m_final.rel is the input for DESeq
#sigASVs_6m is the output from DESeq
#6 MONTHS
sigASVs_6m = sigASVs_6m %>% arrange(log2FoldChange) %>% mutate(Genus = factor(Genus,levels=.$Genus))

DESeq_agp_6m.compositional<- microbiome::transform(agp_infant_6m_final.rel, "compositional")

DESeq_agp_6m_psmelt <- psmelt(DESeq_agp_6m.compositional)

nrow(DESeq_agp_6m_psmelt)
ncol(DESeq_agp_6m_psmelt)
colnames(DESeq_agp_6m_psmelt)
view(DESeq_agp_6m_psmelt)

#renaming OTU column to ASV
colnames(DESeq_agp_6m_psmelt)[1] = "ASV"

DESeq_agp_6m_psmelt_sig <- DESeq_agp_6m_psmelt %>% filter(ASV %in% sigASVs_6m$ASV[sigASVs_6m$padj < 0.05])

nrow(DESeq_agp_6m_psmelt_sig)

#Adding Pseudocount
DESeq_agp_6m_psmelt_sig1 = DESeq_agp_6m_psmelt_sig %>% mutate(Abundance = Abundance + min(.$Abundance[.$Abundance>0]))

DESeq_agp_6m_psmelt_sig1 %>% 
  ggplot(aes(Genus,Abundance,fill=agp_clin)) +
  geom_boxplot() +
  scale_y_log10(expand = expansion(mult = 0.1)) +
  theme_classic(base_size = 16) +
  scale_fill_discrete(labels=c('High', 'Low')) +
  labs(fill='AGP Level') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggpubr::stat_compare_means(aes(group=agp_clin),label='p.signif',size=5)








