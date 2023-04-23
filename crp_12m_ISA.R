#ISA 12 months
library(phyloseq)
library(indicspecies)
library(tidyverse)
#Setting seed
set.seed(1)

#### Make list of datasets to iterate over ####

# Dummy data, with different abundance cutoffs just for fun:
comp = microbiome::transform(crp_infant_12m_final,'compositional')

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#### Create function, set the default value of 'type' as crp ####

# We basically copy your whole code into here, change out crp_infant_12m_final for 'df', and we use 'type' to control parameters that are specific to crp vs CRP analyses
# I've also used ctrl+f to remove 'crp_' from the below code, just so that it's not confusing for others reading the code
run_indicator_analysis = function(df, type = 'crp'){
  
  # Set variable values that you can use to test the function, then comment these lines out once your code's ready to run
  # df = crp_infant_12m_final; type = 'crp'
  
  # infant_12m_final is changed to 'df'
  genus_12m <- tax_glom(df, "Genus", NArm= FALSE)
  genus_RA_12m <- transform_sample_counts(genus_12m, fun=function(x) x/sum(x))
  
  # The first example of using 'type' to define different code for crp vs CRP
  if(type=='crp'){
    isa_12m <- multipatt(t(otu_table(genus_RA_12m)), cluster = sample_data(genus_RA_12m)$crp_median)
  } else if (type=='CRP'){
    isa_12m <- multipatt(t(otu_table(genus_RA_12m)), cluster = sample_data(genus_RA_12m)$crp_clin)
  } else { # Juuust in case
    print('ERROR - unknown "type"')
  }
  
  taxtable_12m_crp <- tax_table(df) %>% as.data.frame() %>% rownames_to_column(var="ASV")
  
  isa_12m$sign = isa_12m$sign %>% rownames_to_column(var="ASV") %>%
    left_join(taxtable_12m_crp) %>% filter(p.value<0.05)
  
  # INDICATOR SPECIES RELATIVE ABUNDANCE
  
  #12 MONTHS
  # Had to change 12m.compositional to compositional.12m since names can't start with numbers
  compositional.12m <- df %>% tax_glom('Genus') %>%  microbiome::transform("compositional")
  
  ISA_12m_psmelt <- psmelt(compositional.12m)
  
  #renaming OTU column to ASV
  colnames(ISA_12m_psmelt)[1] = "ASV"
  
  ISA_12m_psmelt_sig <- ISA_12m_psmelt %>% filter(ASV %in% isa_12m$sign$ASV[isa_12m$sign$p.value<0.05])
  
  #Adding Pseudocount
  ISA_12m_psmelt_sig1 = ISA_12m_psmelt_sig %>% mutate(Abundance = Abundance + min(.$Abundance[.$Abundance>0]))
  
  # Assign plot to 'p' so that we can save it to our output
  # We need 'type' since it differs between analyses
  if(type=='crp'){
    p = ISA_12m_psmelt_sig1 %>% 
      ggplot(aes(Genus,Abundance,fill=crp_median)) +
      geom_boxplot() +
      scale_y_log10(expand = expansion(mult = 0.1)) +
      theme_classic(base_size = 16) +
      scale_fill_discrete(labels=c('High', 'Low')) +
      labs(fill='CRP Level') +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggpubr::stat_compare_means(aes(group=crp_median),label='p.signif',size=6)
  } else if(type=='CRP'){
    p = ISA_12m_psmelt_sig1 %>% 
      ggplot(aes(Genus,Abundance,fill=crp_clin)) +
      geom_boxplot() +
      scale_y_log10(expand = expansion(mult = 0.1)) +
      theme_classic(base_size = 16) +
      scale_fill_discrete(labels=c('High', 'Low')) +
      labs(fill='CRP Level') +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggpubr::stat_compare_means(aes(group=crp_clin),label='p.signif',size=6)
  } else {
    print('Error = unknown "type"')
  }
  
  # All the stuff we want to save in our output
  return(list(ind_output = isa_12m, plot = p))
}

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#### Run the indicator code on each of your datasets of interest ####

df.ind = run_indicator_analysis(crp_infant_12m_final,type='crp')
print(df.ind$plot)
ggsave('ind_crp_12m.jpeg',height=6, width = 5)