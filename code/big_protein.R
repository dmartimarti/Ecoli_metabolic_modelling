library(tidyverse)
library(glue)
library(cowplot)

theme_set(theme_cowplot(14))




gr2199 = gene_PA_long %>% 
  filter(Gene == 'group_2199')

invA = gene_PA_long %>% 
  filter(Gene == 'mxiA~~~invA~~~flhA_3~~~flhA_4~~~flhA_2~~~flhA_1')


gr2199 %>% 
  rename(gr2199  =  presence) %>% 
  select(Genome, gr2199) %>% 
  left_join(
    invA %>% 
      rename(inv = presence) %>% 
      select(Genome, inv)
  ) %>% 
  mutate(total_pres = gr2199 + inv) %>% 
  filter(total_pres != 0) %>% 
  mutate(total_pres = factor(total_pres)) %>% 
  filter(inv != 0) %>% 
  count(total_pres) %>% 
  mutate(total_pres  = case_when(total_pres == 1 ~ 'Only one gene in genome',
                                 TRUE ~ 'Both genes in genome')) %>% 
  rename(Group=total_pres) %>% 
  ggplot(aes(y = Group, x = n)) +
  geom_segment(aes(y=Group, yend=Group, x=0, xend=n)) +
  geom_point(shape = 21, color = 'red', fill = alpha('orange', 0.8), size = 5)
  
  

# get a list of strains that contain the gene and filter metadata with them

metadata %>% filter( Genome %in%  
                       (gr2199 %>% filter(presence == 1) %>% 
                          pull(Genome))) %>% 
  write_csv('prot_structures_esm/big_protein/strains_with_group2199.csv')


genomes_gr2199 = gr2199 %>% filter(presence == 1) %>% 
  pull(Genome)
genomes_gr2199_ID = metadata %>% filter(Genome %in% genomes_gr2199) %>% 
  pull(ID)

# check wether the proportion commensal/pathogen is different with and without
# the gene
# it isn't!

metadata %>% 
  filter(Discard == 'No') %>% 
  mutate(big_gene = case_when(Genome %in% genomes_gr2199 ~ 'Yes',
                              TRUE ~ 'No')) %>% 
  drop_na(Broadphenotype) %>% 
  filter(Broadphenotype != "Laboratory strain") %>% 
  filter(Broadphenotype != 'Unknown') %>% 
  distinct(Genome, .keep_all = T) %>% 
  count(Broadphenotype, big_gene) %>% 
  pivot_wider(names_from = big_gene, values_from = n) %>% 
  select(-Broadphenotype) %>% 
  chisq.test()
  
  
metadata %>% 
  filter(Discard == 'No') %>% 
  mutate(big_gene = case_when(Genome %in% genomes_gr2199 ~ 'Yes',
                              TRUE ~ 'No')) %>% 
  drop_na(Broadphenotype) %>% 
  filter(Broadphenotype != "Laboratory strain") %>% 
  filter(Broadphenotype != 'Unknown') %>% 
  distinct(Genome, .keep_all = T) %>% 
  count(Broadphenotype, big_gene) %>% 
  group_by(big_gene) %>% 
  mutate(total = sum(n),
         prop = n / total) %>% 
  ggplot(aes(x = Broadphenotype, y = n, fill = big_gene)) +
  geom_bar(position="fill", stat="identity")




## difference in the worm brightness?

library(readxl)
ALL_worm_FC <- read_excel("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/worm_imaging/ALL_worm_FC.xlsx")

worm_gr2199 = ALL_worm_FC %>% 
  filter(ID %in% (metadata %>% 
           filter(Discard == 'No') %>% 
           pull(ID))) %>% 
  distinct(ID, .keep_all = T) %>% 
  mutate(big_gene = case_when(ID %in% genomes_gr2199_ID ~ 'Yes',
                              TRUE ~ 'No')) 

worm_gr2199 %>% 
  ggplot(aes(x = big_gene, y = Mean_FC, fill = big_gene)) +
  geom_boxplot(show.legend = F) +
  geom_point(position = position_jitterdodge(), show.legend = F)
  

worm_gr2199 %>% 
  rstatix::t_test(Mean_FC ~ big_gene, detailed = T)







