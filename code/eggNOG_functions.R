# libraries

library(tidyverse)
library(readr)
library(tidymodels)
library(viridis)
library(ggpubr)


# read data ---------------------------------------------------------------


# read the eggNOG annotations

eggnog =  read_delim("tables/eggnog_PG_annotation.tsv", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

eggnog = eggnog %>% 
  rename(Gene = `#query`)

eggnog %>% 
  write_csv('tables/eggnog_PG_annotation.csv')

# load the full list of genes and PA in our genomes
gene_PA = read_delim("tables/gene_presence_absence.Rtab", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

gene_PA = gene_PA %>% 
  mutate(sum = rowSums(select(., -Gene)),
         sum = sum/max(sum),
         group = case_when(sum > 0.99 ~ 'core',
                           sum <= 0.99 & sum > 0.95 ~ 'soft_core',
                           sum <= 0.95 & sum > 0.15 ~ 'shell',
                           sum <= 0.15 ~ 'cloud'), 
         .before = Gene)



## add info for the COG categories
# Create the first and second column vectors
COG_category = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", 
                 "K", "L", "M", "N", "O", "P", "Q", "T", "U", 
                 "Y", "Z", "R", "S", "V")

COG_description = c("RNA processing and modification", 
                    "Chromatin Structure and dynamics", 
                    "Energy production and conversion", 
                    "Cell cycle control and mitosis", 
                    "Amino Acid metabolis and transport",
                    "Nucleotide metabolism and transport", 
                    "Carbohydrate metabolism and transport", 
                    "Coenzyme metabolis", "Lipid metabolism", 
                    "Tranlsation", "Transcription", "Replication and repair", 
                    "Cell wall/membrane/envelop biogenesis", "Cell motility", 
                    "Post-translational modification, protein turnover, chaperone functions", 
                    "Inorganic ion transport and metabolism", "Secondary Structure", 
                    "Signal Transduction", "Intracellular trafficing and secretion", 
                    "Nuclear structure", "Cytoskeleton", 
                    "General Functional Prediction only", "Function Unknown",
                    "Defense mechanism")


# Create a data frame with the first and second column vectors
COG_descr_df = tibble(COG_category, COG_description)



# merge datasets ----------------------------------------------------------


# join the table and see how many gene/proteins were annotated somehow with eggnog
eggnog_PA = eggnog %>% 
  select(Gene, Preferred_name, COG_category, GOs, KEGG_ko, KEGG_Pathway, 
         KEGG_Module, KEGG_Reaction) %>% 
  # separate COG categories
  separate_rows(COG_category, sep = '') %>%
  filter(COG_category != '') %>% # filter crap away
  # add the COG info from COG categories
  left_join(COG_descr_df) %>% 
  group_by(Gene) %>% 
  mutate(COG_category =  paste(COG_category, collapse = ','),
         COG_description = paste(COG_description, collapse = ',')) %>% 
  ungroup() %>% 
  distinct(Gene, .keep_all = T) %>% 
  select(Gene, Preferred_name, COG_category, COG_description, 
         everything()) %>% 
  mutate(annotated = 'Yes') %>% 
  full_join(gene_PA %>% 
              select(Gene, group)) %>% 
  mutate(annotated = case_when(annotated == 'Yes' ~ 'Yes', 
                               TRUE ~ 'No')) %>% 
  select(Gene, group, annotated, COG_category, COG_description,
         everything())

eggnog_PA = eggnog_PA %>% 
  mutate(COG_category = case_when(COG_category == '-' ~ 'Not annotated',
                                  is.na(COG_category) ~ 'Not annotated',
                                  TRUE ~ COG_category))

eggnog_PA = eggnog_PA %>% 
  mutate(Preferred_name = case_when(is.na(Preferred_name) ~ '-',
                                    TRUE ~ Preferred_name)) 

eggnog_PA %>% 
  write_csv("tables/eggnog_annotation_fixed.csv")

# genome groups comparison ------------------------------------------------


# plot the proportion of genes annotated in eggnog
eggnog_PA %>% 
  group_by(group, annotated) %>% 
  count() %>% 
  drop_na(group) %>% 
  group_by(group) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup %>% 
  mutate(group = factor(group, levels = c('core', 'soft_core', 
                                          'shell', 'cloud'))) %>% 
  ggplot(aes(x = group, y = prop, fill = annotated)) +
  geom_col(position="fill") +
  scale_fill_viridis(discrete = T, option = "E", direction = -1) +
  labs(x = 'Genome group',
       y = 'Proportion',
       fill='Annotated\nin eggNOG?', 
       caption = 
         'Proportion of gene families annotated with eggNOG')

ggsave("exploration/eggnog/eggnog_annot_props.pdf", 
       height = 6, width = 6)


# proportion of known and unknown gene/protein functions
eggnog_PA %>% 
  separate_rows(COG_category, sep = ',') %>% 
  mutate(known = case_when(
    COG_category %in% c('S', 'Not annotated') ~ 'Uknown',
    TRUE ~ 'Known'
  )) %>% 
  group_by(group) %>% 
  count(known) %>% 
  group_by(group) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup %>% 
  mutate(group = factor(group, levels = c('core', 'soft_core', 
                                          'shell', 'cloud'))) %>% 
  ggplot(aes(x = group, y = prop, fill = known)) +
  geom_col(position="fill") +
  scale_fill_viridis(discrete = T, option = "E") + 
  labs(x = 'Genome group',
       y = 'Proportion',
       fill='Annotation', 
       caption = 
         'Gene families with known COG annotation')

ggsave("exploration/eggnog/eggnog_known_functions.pdf", 
       height = 6, width = 6)


# COG proportions
eggnog_PA %>% 
  separate_rows(COG_category, sep = ',') %>% 
  mutate(known = case_when(
    COG_category %in% c('S', 'Not annotated') ~ 'Uknown',
    TRUE ~ 'Known'
  )) %>% 
  filter(known == 'Known') %>% 
  group_by(group, COG_category) %>% 
  count() %>% 
  group_by(group) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup %>% 
  mutate(group = factor(group, levels = c('core', 'soft_core', 
                                          'shell', 'cloud'))) %>% 
  ggplot(aes(x = prop, y = fct_reorder(COG_category, prop), fill = COG_category)) +
  geom_col(position="dodge", show.legend = FALSE) +
  facet_wrap(~group) +
  # scale_fill_viridis(discrete = T, option = "E") +
  labs(x = 'Proportion',
       y = 'COG category',
       fill='COG category')


ggsave("exploration/eggnog/eggnog_COG_categories.pdf", 
       height = 8, width = 9)


# heatmap
library(ComplexHeatmap)

eggnog_mat = eggnog_PA %>% 
  separate_rows(COG_category, sep = ',') %>% 
  mutate(known = case_when(
    COG_category %in% c('S', 'Not annotated') ~ 'Uknown',
    TRUE ~ 'Known'
  )) %>% 
  filter(known == 'Known') %>% 
  group_by(group, COG_category) %>% 
  count() %>% 
  group_by(group) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup %>% 
  select(group, COG_category, prop) %>% 
  pivot_wider(names_from = group, values_from = prop) %>% 
  replace_na(list(cloud = 0, core = 0, shell = 0, soft_core = 0)) %>% 
  column_to_rownames('COG_category') %>% 
  as.matrix()

max_val = max(eggnog_mat)

chisq.test(eggnog_mat)

# Calculate the IQR for each row of the matrix
iqr = apply(eggnog_mat, MARGIN = 1, function(x) {
  IQR(x)
})

# Calculate the quartile values for each row of the matrix
quartiles = apply(eggnog_mat, MARGIN = 1, function(x) {
  quantile(x, probs = c(0.25, 0.5, 0.75))
})


# calculating enrichment of categories if they are larger than
# (IQR*1.5) + upper_quartile
enrich_mat = rbind(iqr, quartiles) %>% 
  t %>% 
  as_tibble() %>% rownames_to_column('category') %>% 
  rename(lower = `25%`,
         median = `50%`,
         upper = `75%`) %>% 
  mutate(category = colnames(quartiles),
         upp_limit = (1.5 * iqr) + upper) %>% 
  bind_cols(eggnog_mat) %>% 
  mutate(core = case_when(core > upp_limit ~ '*',
                          TRUE ~ ''),
         cloud = case_when(cloud > upp_limit ~ '*',
                           TRUE ~ ''),
         shell = case_when(shell > upp_limit ~ '*',
                          TRUE ~ ''),
         soft_core = case_when(soft_core > upp_limit ~ '*',
                           TRUE ~ '')) %>% 
  select(category, cloud:soft_core) %>% 
  column_to_rownames('category') %>% 
  as.matrix()



ht = eggnog_mat %>% 
  Heatmap(
    col = colorRamp2(c(0, max_val * 1.05), c("white", "darkblue")),
    name = 'Proportion',
    column_names_gp = gpar(fontsize = 20),
    column_names_rot = 0,
    cell_fun = function(j,i,x,y,width,heigh,fill) {
      grid.text(sprintf("%.1s", enrich_mat[i,j]),
                x, y, gp = gpar(fontsize = 10))
    }
  )  

ht


# boxplot of proportions
eggnog_PA %>% 
  separate_rows(COG_category, sep = ',') %>% 
  mutate(known = case_when(
    COG_category %in% c('S', 'Not annotated') ~ 'Uknown',
    TRUE ~ 'Known'
  )) %>% 
  filter(known == 'Known') %>% 
  group_by(group, COG_category) %>% 
  count() %>% 
  group_by(group) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup %>% 
  select(group, COG_category, prop) %>% 
  ggplot(aes(x=COG_category, y = prop, fill = COG_category)) +
  geom_boxplot(show.legend = F) +
  geom_point(aes(fill = group), shape = 21, size = 4, color = 'black',
             position = position_jitterdodge(dodge.width=0.5))


ggsave("exploration/eggnog/COG_categories_boxplot_prop.pdf", 
       height = 6, width = 7)



# phylogroups core genome -------------------------------------------------



gene_PA_long = gene_PA %>% 
  pivot_longer(-sum:-Gene, 
               names_to = 'Genome', 
               values_to = 'presence')

unique(gene_PA_long$Genome)

metadata

gene_PA_metadata = gene_PA_long %>% 
  select(-sum) %>% 
  left_join(metadata %>% 
              filter(Origin %in% c('ECOREF', 'AUS'),
                     Discard == 'No') %>% 
              select(Genome, ID, Origin, phylogroup, Broadphenotype)) %>% 
  drop_na(phylogroup ) %>% 
  select(-Origin) %>% 
  mutate(phylogroup = case_when(phylogroup == 'E or cladeI' ~ 'E',
                                TRUE ~ phylogroup)) %>% 
  filter(phylogroup != 'G')




## phylogroup A ------------------------------------------------------------


phylo_core = function(phylo = 'A') {
  # total members in the phylogroup 
  n_A = gene_PA_metadata %>% 
    filter(phylogroup == phylo) %>% 
    distinct(Genome) %>% count() %>% pull(n)
  
  size = gene_PA_metadata %>% 
    filter(phylogroup == phylo) %>% 
    group_by(Gene) %>% 
    summarise(pres_A = sum(presence),
              sum = pres_A / n_A) %>% 
    filter(pres_A != 0) %>% 
    mutate(group = case_when(sum > 0.99 ~ 'core',
                             sum <= 0.99 & sum > 0.95 ~ 'soft_core',
                             sum <= 0.95 & sum > 0.15 ~ 'shell',
                             sum <= 0.15 ~ 'cloud')) %>% 
    group_by(group) %>% 
    count() %>% 
    ungroup
  
  return(size)
}

A_size = phylo_core(phylo = 'A') %>% rename(A = n) 
B1_size = phylo_core(phylo = 'B1') %>% rename(B1 = n)
B2_size = phylo_core(phylo = 'B2') %>% rename(B2 = n) 
C_size = phylo_core(phylo = 'C') %>% rename(C = n)
D_size = phylo_core(phylo = 'D') %>% rename(D = n) 
E_size = phylo_core(phylo = 'E') %>% rename(E = n)
F_size = phylo_core(phylo = 'F') %>% rename(F = n) 

pangenome_sizes = A_size %>% 
  left_join(B1_size) %>% 
  left_join(B2_size) %>% 
  left_join(C_size) %>% 
  left_join(D_size) %>% 
  left_join(E_size) %>% 
  left_join(F_size) %>% 
  replace_na(list(C = 0))

# include the total gene count
pangenome_sizes = pangenome_sizes %>% 
  summarise(across(A:F, sum)) %>% 
  mutate(group = 'total_genes', .before='A') %>% 
  bind_rows(pangenome_sizes)


phylos = unique(colnames(pangenome_sizes)[-1])
member_size = c()
for (phylo in phylos) {
  n = gene_PA_metadata %>% 
    filter(phylogroup == phylo) %>% 
    distinct(Genome) %>% count() %>% pull(n)
  member_size = c(member_size, n)
}

member_size
# add the info of the member size for each phylogroup
pangenome_sizes = pangenome_sizes %>% 
  rbind(c('members', member_size))


pangenome_sizes %>% 
  pivot_longer(
    -group,
    names_to = 'phylogroup',
    values_to = 'genes') %>% 
  mutate(genes = as.integer(genes)) %>% 
  filter(!(group %in% c('members','total_genes'))) %>% 
  ggplot(aes(x = phylogroup, y = genes, fill = group)) +
  geom_col(position = 'dodge')

ggsave("exploration/eggnog/phylogroups_genes_cats.pdf", 
       height = 6, width = 7)



# correlations between pair of values -------------------------------


pan_corr = pangenome_sizes %>% 
  pivot_longer(
    -group,
    names_to = 'phylogroup',
    values_to = 'genes') %>% 
  pivot_wider(names_from = group, values_from = genes) %>% 
  mutate(across(total_genes:members, as.integer))

# members vs core
pan_corr %>% 
  ggplot(aes(x = core, y = members)) +
  geom_smooth(method = 'lm') + 
  geom_point() +
  stat_cor(method = "pearson", label.x = 3000, label.y = 310)


ggsave("exploration/eggnog/correlations/core_members_pearson.pdf", 
       height = 6, width = 7)

# members vs soft
pan_corr %>% 
  ggplot(aes(x = total_genes, y = cloud)) +
  geom_smooth(method = 'lm') + 
  geom_point() +
  stat_cor(method = "pearson")


ggsave("exploration/eggnog/correlations/cloud_members_pearson.pdf", 
       height = 6, width = 7)




pdf("exploration/eggnog/correlations/pairwise_pearson.pdf",
    height = 6, width = 6.5)
PerformanceAnalytics::chart.Correlation(pan_corr[,2:7], histogram=TRUE, pch=19)
dev.off()





# alleles present in pangenome --------------------------------------------


eggnog_PA %>% 
  filter(str_detect(Gene, 'group_')) %>% 
  filter(Preferred_name != '-') %>% view



eggnog_PA %>% 
  filter(!(str_detect(Gene, 'group_'))) %>%
  mutate(allele = case_when(str_detect(Gene, '~~~') ~ 'allele',
                            TRUE ~ 'not_allele')) %>% 
  filter()
  separate_rows(Gene, sep = '~~~') %>% 
  mutate(Gene_unique = str_replace(Gene, '_[:digit:]', ""),
         Gene_unique = str_replace(Gene_unique, "[:digit:]", ''), 
         .before = 'group') 


eggnog_PA %>% 
  filter(!(str_detect(Gene, 'group_'))) %>%
  mutate(allele = case_when(str_detect(Gene, '~~~') ~ 'allele',
                            TRUE ~ 'not_allele'))


























