
# libraries ---------------------------------------------------------------



library(tidyverse)
library(readr)
library(cowplot)
library(here)
library(ComplexHeatmap)
library(circlize)
library(broom)
library(ggsankey)
# install.packages("vcd")
library(vcd)
library(glue)
library(openxlsx)


theme_set(theme_cowplot(15))



# load data ---------------------------------------------------------------



# get the genome names
gnm_names = list.files(path = ".", pattern = ".tbl") %>% 
  str_sub(start = 1, end = -18)

# get the files names
files_list = list.files(path = ".", pattern = ".tbl")


# loop over the files and get them in a huge table
genome_paths = tibble()
for (i in 1:length(gnm_names)){
  temp = read_delim(files_list[i], 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE, skip = 2) %>% 
    mutate(Genome = gnm_names[i])
  
  genome_paths = bind_rows(genome_paths, temp)
}

genome_paths %>% 
  distinct(ID, .keep_all = TRUE) %>% view


genome_paths %>% 
  filter(ID == '|12DICHLORETHDEG-PWY|') %>% View

# Metadatas ####
## genome metadata  ####

library(readxl)
metadata = read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx") %>% 
  mutate(Genome = str_sub(fasta, end = -7), .before = Strainname) %>% 
  filter(Genome %in% gnm_names) %>% 
  filter(phylogroup != 'cladeI') %>% 
  distinct(Genome, .keep_all = TRUE)

## pathways-reactions metadata ####

metacyc_paths = read_delim("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/gapseq/All-reactions-of-MetaCyc_extended.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) %>% 
  rename(ec_number = `EC-Number`,
         ID = `Object ID`,
         reaction = Reaction,
         substrates = `Substrates`,
         common_name = `Common-Name`,
         names = Names,
         spontaneous = `Spontaneous?`,
         systematic_name= `Systematic-Name`,
         pathway = `In-Pathway`)


pan_reactions = genome_paths %>% 
  arrange(desc(Completeness)) %>% 
  distinct(ID, .keep_all = T) %>% 
  separate_rows(ReactionsFound, sep = ' ') %>% 
  distinct(ReactionsFound) %>% 
  pull(ReactionsFound)

metacyc_paths_pangenome = metacyc_paths %>% 
  filter(reaction %in% pan_reactions)  %>% 
  separate_rows(substrates, 
                sep = ' // ') %>% 
  rename(met_name = substrates)


## ModelSeed comps ####

modelSeed_reactions = read_delim("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/reactions.tsv", 
           delim = "\t", escape_double = FALSE, 
           trim_ws = TRUE)


modelSeed_comps = read_delim("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)



# prepare the reactions dataset
modelSeed_reac_pan = modelSeed_reactions %>% 
  separate_rows(aliases, sep = '\\|') %>% 
  filter(str_detect(aliases, 'MetaCyc:')) %>% 
  mutate(aliases = str_sub(aliases, start = 10)) %>% 
  filter(aliases %in% pan_reactions) %>% 
  select(-notes, -source) 


modelSeed_reac_pan %>% 
  separate_rows(compound_ids, sep = ';')




# how many genomes are we lacking ? 

# filter metadata for the bugs we want to analyse: AUS and ECOREF
meta_filt = metadata %>% 
  mutate(Genome = str_sub(fasta, start = 1, end = -7), .before = ID) %>% 
  filter(Discard == 'No') %>% 
  filter(Origin %in% c('AUS', 'ECOREF')) %>% 
  distinct(Genome, .keep_all = T) 


# save the missing genome list in a csv file
meta_filt %>% 
  filter(!(Genome %in% gnm_names)) %>% 
  write_csv('../missing_genomes.csv')
  





# exploration plots -------------------------------------------------------

# first filter the genome list we want to analyse

genome_paths = genome_paths %>% filter(Genome %in% meta_filt$Genome)


# calculate how many pathways per genome are in our dataset


## histograms ####

genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  group_by(Genome) %>% 
  count() %>% 
  filter(n > 400) %>% 
  ggplot(aes(n)) +
  geom_histogram(color = 'black', fill = 'grey50', bins = 30) +
  labs(
    x = 'Number of complete pathways in genome',
    y = 'Number of genomes'
  )

ggsave("../exploration/number_of_paths.pdf", height = 7, width = 9)


genome_paths %>% 
  filter(Completeness > 0) %>% 
  select(Name, Genome) %>% 
  group_by(Genome) %>% 
  count() %>% 
  filter(n > 400) %>%
  ggplot(aes(n)) +
  geom_histogram(color = 'black', fill = 'grey50', bins = 30) +
  labs(
    x = 'Number of pathways in genome',
    y = 'Number of genomes',
    caption = 'Pathways not complete (<100% of coverage) are also included in this plot'
  )
ggsave("../exploration/number_of_paths_all.pdf", height = 7, width = 9)

## paths per phylogroup ####
# distinct paths per phylogroup

genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  left_join(meta_filt %>% 
              select(Genome, phylogroup)) %>% 
  group_by(Genome, phylogroup) %>% 
  count() %>% 
  filter(n > 400) %>%
  ggplot(aes(n, fill = phylogroup
             )) +
  geom_histogram(color = 'black', bins = 30) +
  labs(
    x = 'Number of pathways in genome',
    y = 'Number of genomes'
  )
 
ggsave("../exploration/number_of_paths_phylogroups.pdf", height = 7, width = 9)

genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  left_join(meta_filt %>% 
              select(Genome, phylogroup)) %>% 
  group_by(Genome, phylogroup) %>% 
  filter(phylogroup != 'cladeI') %>% 
  count() %>% 
  filter(n > 400) %>%
  ggplot(aes(n, fill = phylogroup
  )) +
  geom_density(color = 'black', alpha = 0.9) +
  labs(
    x = 'Number of pathways in genome',
    y = 'Number of genomes'
  ) + 
  facet_wrap(~phylogroup, ncol = 1)

ggsave("../exploration/paths_density_phylogroups.pdf", height = 12, width = 6)

# pathways heatmaps ####
## complete paths matrix ####
paths_pa = genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0)

pahts_names = names(paths_pa)
gnm_names = paths_pa$Genome

paths_matrix = paths_pa %>% 
  select(-Genome) %>% 
  as.matrix() 
  
rownames(paths_matrix) = gnm_names


col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))
col_fun(seq(-3, 3))
Heatmap(paths_matrix, 
        col = col_fun,
        name = "Completeness")


quartz.save(file = '../exploration/Pathways_heatmap.pdf',
            type = 'pdf', height = 70, width = 90)


Heatmap(paths_matrix, 
        col = col_fun,
        name = "Completeness",
        show_row_names = FALSE,
        show_column_names = FALSE)

quartz.save(file = '../exploration/Pathways_heatmap_noNames.pdf',
            type = 'pdf', height = 9, width = 12)


## difference paths matrix  ####

paths_pa = genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0)

pahts_names = names(paths_pa)
gnm_names = paths_pa$Genome

paths_matrix = paths_pa %>% 
  select(-Genome) %>% 
  as.matrix() 

rownames(paths_matrix) = gnm_names

# total genomes
tot_gnms = dim(paths_matrix)[1]

# filter the paths that are NOT 100% complete in at least 99% of genomes
# (core pathway)
redux_matrix = paths_matrix[,colSums(paths_matrix) < tot_gnms*0.99]
dim(redux_matrix)

### differnetial pathways ####
diff_pathways = colnames(redux_matrix)

col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))
Heatmap(redux_matrix, 
        col = col_fun,
        name = "Completeness")

quartz.save(file = '../exploration/Pathways_heatmap_differential.pdf',
            type = 'pdf', height = 50, width = 70)

Heatmap(redux_matrix, 
        col = col_fun,
        name = "Completeness",
        show_row_names = FALSE,
        show_column_names = FALSE)

quartz.save(file = '../exploration/Pathways_heatmap_differential_noNames.pdf',
            type = 'pdf', height = 6, width = 9)



## all paths matrix ####

paths_pa = genome_paths %>% 
  filter(Completeness > 0) %>% 
  mutate(presence = Completeness/100) %>% 
  select(Name, Genome, presence) %>% 
  # mutate(presence = 1) %>% 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0)

pahts_names = names(paths_pa)
gnm_names = paths_pa$Genome

paths_matrix = paths_pa %>% 
  select(-Genome) %>% 
  as.matrix() 

rownames(paths_matrix) = gnm_names

Heatmap(paths_matrix, 
        col = col_fun,
        name = "Completeness")

quartz.save(file = '../exploration/ALL_Pathways_heatmap.pdf',
            type = 'pdf', height = 70, width = 190)

Heatmap(paths_matrix, 
        col = col_fun,
        name = "Completeness",
        show_row_names = FALSE,
        show_column_names = FALSE)

quartz.save(file = '../exploration/ALL_Pathways_heatmap_noNames.pdf',
            type = 'pdf', height = 9, width = 12)






# metabolites -------------------------------------------------------

# read metabolites extracted from the R objects (see the other script)
metab = read_csv("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/gapseq/model_metabolites.csv")


gnm_valid = meta_filt %>% pull(Genome)

unique(metab$comp)

gnm_remove = c('SPC_3.1','OP50',
               'NT12226_337','NT12229','NT12169_305',
               '157')

metab %>% 
  filter(!(model %in% gnm_remove)) %>% 
  filter(model %in% gnm_valid) %>% 
  group_by(model) %>% 
  count() %>% 
  ggplot(aes(n)) +
  geom_histogram(color = 'black', fill = 'grey50', bins = 30) +
  labs(
    x = 'Molecules',
    y = 'Count'
  )

ggsave("../exploration/molecules_histogram.pdf", height = 7, width = 9)
  



metab %>% distinct(comp)

## periplasm metabs ####

# simple version with ggplot
metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'p0') %>% 
  mutate(val = 1) %>% 
  ggplot(aes(met_name, model, fill = val)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))


# generate a table with the differntial metabolites
metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'p0') %>% 
  select(met_name, model) %>% 
  mutate(presence = 1) %>% 
  count(met_name) %>% 
  arrange(n) %>% 
  # filter(n < ((736 * 0.1))) %>% 
  write_csv('../exploration/tables/metabs_periplasm.csv')

metacyc_paths_pangenome

# complexHeatmap
mets_pa = metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'p0') %>% 
  select(met_name, model) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = met_name, values_from = presence, values_fill = 0)

mets_names = names(mets_pa)
gnm_names = mets_pa$model

mets_matrix = mets_pa %>% 
  select(-model) %>% 
  as.matrix() 

rownames(mets_matrix) = gnm_names


col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))

Heatmap(mets_matrix, 
        col = col_fun,
        name = "Presence")


quartz.save(file = '../exploration/Metabs_periplasm_heatmap.pdf',
            type = 'pdf', height = 80, width = 9)





## external metabs ####


metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'e0') %>% 
  mutate(val = 1) %>% 
  ggplot(aes(met_name, model, fill = val)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))

# generate a table with the differntial metabolites
metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'e0') %>%
  group_by(met_id) %>% 
  count(met_name) %>% 
  arrange(n) %>% 
  write_csv('../exploration/tables/metabs_ext.csv')

metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'e0') %>% 
  select(met_name, model) %>% 
  mutate(presence = 1) %>% 
  count(met_name) %>% 
  arrange(n) %>% 
  full_join(metacyc_paths_pangenome) %>% 
  write_csv('../exploration/tables/metabs_ext_descriptions.csv')




# complexHeatmap
mets_pa = metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'e0') %>% 
  select(met_name, model) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = met_name, values_from = presence,
              values_fn = mean, values_fill = 0)

mets_names = names(mets_pa)
gnm_names = mets_pa$model

mets_matrix = mets_pa %>% 
  select(-model) %>% 
  as.matrix() 

rownames(mets_matrix) = gnm_names

### reduce the matrix ####
# total genomes
tot_gnms = dim(mets_matrix)[1]

# filter the paths that are NOT 100% complete in at least 99% of genomes
# (core pathway)
redux_matrix = mets_matrix[,colSums(mets_matrix) < tot_gnms*0.99]
dim(redux_matrix)

# save the variables for the PCA later
metab_ext_matrix = mets_matrix
metab_ext_redux_matrix = redux_matrix


col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))

Heatmap(metab_ext_matrix, 
        col = col_fun,
        name = "Presence")



quartz.save(file = '../exploration/Metabs_ext_heatmap.pdf',
            type = 'pdf', height = 80, width = 30)



Heatmap(redux_matrix, 
        col = col_fun,
        name = "Presence",
        show_row_names = F)


quartz.save(file = '../exploration/Metabs_ext_redux_nonames_heatmap.pdf',
            type = 'pdf', height = 7, width = 8)


## cytoplasm metabs ####

metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'c0') %>% 
  mutate(val = 1) %>% 
  ggplot(aes(met_name, model, fill = val)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))


# generate a table with the differential metabolites
metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'c0') %>%
  group_by(met_id) %>% 
  count(met_name) %>% 
  ungroup %>% 
  arrange(n) %>% 
  write_csv('../exploration/tables/metabs_cyt.csv')


metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'c0') %>% 
  select(met_name, model) %>% 
  mutate(presence = 1) %>% 
  count(met_name) %>% 
  arrange(n) %>% 
  full_join(metacyc_paths_pangenome) %>% 
  # filter(n < ((736 * 0.1))) %>% 
  write_csv('../exploration/tables/metabs_cyt_description.csv')

  

# complexHeatmap
mets_pa = metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'c0') %>% 
  select(met_name, model) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = met_name, values_from = presence,
              values_fn = mean, values_fill = 0)

mets_names = names(mets_pa)
gnm_names = mets_pa$model

mets_matrix = mets_pa %>% 
  select(-model) %>% 
  as.matrix() 

rownames(mets_matrix) = gnm_names

### reduce the matrix ####
# total genomes
tot_gnms = dim(mets_matrix)[1]

# filter the paths that are NOT 100% complete in at least 99% of genomes
# (core pathway)
redux_matrix = mets_matrix[,colSums(mets_matrix) < tot_gnms*0.99]
dim(redux_matrix)

# save the variables for the PCA later
metab_cyto_matrix = mets_matrix
metab_cyto_redux_matrix = redux_matrix

col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))

Heatmap(redux_matrix, 
        col = col_fun,
        name = "Presence")


quartz.save(file = '../exploration/Metabs_cytop_heatmap.pdf',
            type = 'pdf', height = 80, width = 100)



# PCA of metabolites ------------------------------------------------------

library(broom)


# External mets
# metab_ext_matrix
# metab_ext_redux_matrix

## ext metabolites ####

pca_fit = metab_cyto_matrix %>% 
  as.data.frame() %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale = F)

pca_fit %>%
  augment(redux_matrix) %>% # add original dataset back in
  rename(Genome = `.rownames`) %>% 
  select(Genome,`.fittedPC1`:`.fittedPC5`) %>% 
  left_join(meta_filt) %>% 
  filter(!(phylogroup %in% c('E or cladeI', 'cladeI'))) %>% 
  # filter(phylogroup != 'B2') %>% 
  ggplot(aes(.fittedPC1, .fittedPC3, 
             color = phylogroup, fill = phylogroup)) + 
  geom_point(size = 1.5) +
  # ylim(-2,2) +
  theme_half_open(12) + 
  labs(
    x = 'PC1',
    y = 'PC3',
    caption = ''
    ) +
  stat_ellipse(level=0.95, geom = 'polygon', alpha = 0.3) +
  background_grid()

ggsave("../exploration/PCA_ext_mets_PC1_PC3.pdf", height = 7, width = 9)


## cyto metabolites ####

# cytoplasm mets
# metab_cyto_matrix
# metab_cyto_redux_matrix


pca_fit = metab_cyto_redux_matrix %>% 
  as.data.frame() %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale = F)

pca_fit %>%
  augment(redux_matrix) %>% # add original dataset back in
  rename(Genome = `.rownames`) %>% 
  select(Genome,`.fittedPC1`:`.fittedPC5`) %>% 
  left_join(meta_filt) %>% 
  filter(!(phylogroup %in% c('E or cladeI', 'cladeI'))) %>% 
  # filter(phylogroup != 'B2') %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, 
             color = phylogroup, fill = phylogroup)) + 
  geom_point(size = 1.5) +
  # ylim(-2,2) +
  theme_half_open(12) + 
  labs(
    x = 'PC1',
    y = 'PC2',
    caption = ''
  ) +
  stat_ellipse(level=0.95, geom = 'polygon', alpha = 0.3) +
  background_grid()

ggsave("../exploration/PCA_cyto_mets.pdf", height = 7, width = 9)







# Chi-square test ---------------------------------------------------------

# prepare the matrix
chi_sq_matrix = genome_paths %>% 
  filter(Completeness == 100) %>% 
  # filter(Name %in% diff_pathways) %>%  ## filter the pathways that seem to be different
  left_join(metadata %>% 
              select(Genome, phylogroup)) %>%
  drop_na(phylogroup) %>% 
  mutate(phylogroup = case_when(phylogroup == 'E or cladeI' ~ 'E',
                                TRUE ~ phylogroup)) %>% 
  group_by(Name, phylogroup) %>% 
  count() %>% 
  pivot_wider(names_from = phylogroup, values_from = n, values_fill = 0)


pahts_names = chi_sq_matrix$Name

chi_sq_matrix = chi_sq_matrix %>% 
  ungroup %>% 
  select(-Name) %>% 
  as.data.frame() 

rownames(chi_sq_matrix) = pahts_names


# convert to table
dt = as.table(as.matrix(chi_sq_matrix))



## test table with 5 pathways ####

test_paths = c('ATP biosynthesis', 'isoniazid activation', # present in all
               'Agmatine Transport', 'mannitol cycle', # present in all
               'homotaurine degradation', 'acetoin degradation', # present in all
               'curcumin degradation', 'maltose degradation', # present in a few
               'salmochelin degradation', 'salmochelin biosynthesis') # present in a few

test_dt = dt[test_paths,]

assoc(test_dt, shade = TRUE, las=1,
      labeling_args = list(
        offset_labels = c(left = -4),
        offset_varnames = c(left = 0),
        rot_labels = c(left = 0)))

quartz.save(file = '../exploration/chisq_plots/assoc_test.pdf',
            type = 'pdf', height = 10, width = 7)


mosaicplot(test_dt , shade = TRUE, las=2,
           cex.axis = 01,
           main = "sdlfakjsd")


chisq = chisq.test(test_dt)

chisq$observed

round(chisq$expected,2)
round(chisq$residuals,3)

corrplot::corrplot(chisq$residuals, is.cor = FALSE)


# contribution 
contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)
corrplot::corrplot(contrib, is.cor = FALSE)

## complete dataset

mosaicplot(t(head(dt,20)), shade = TRUE, las=2,
           cex.axis = 01,
           main = "sdlfakjsd")

quartz.save(file = '../exploration/chisq_plots/mosaicplot_diff_paths.pdf',
            type = 'pdf', height = 18, width = 16)


# plot just a subset of the table
assoc(dt, shade = TRUE, 
      las=1,
      labeling_args = list(
        offset_labels = c(left = -4),
        offset_varnames = c(left = 0),
        rot_labels = c(left = 0)))

quartz.save(file = '../exploration/chisq_plots/assoc_all.pdf',
            type = 'pdf', height = 120, width = 12)

## chi sq test ####
chisq = chisq.test(dt)

chisq$observed

round(chisq$expected,2)
residuals = round(chisq$residuals,3)





corrplot::corrplot(residuals[rowSums(abs(residuals)) >15,], 
                   is.cor = FALSE,
                   tl.cex = 0.3)

quartz.save(file = '../exploration/chisq_plots/residuals_15.pdf',
            type = 'pdf', height = 10, width = 5)

# contribution 
contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)
corrplot::corrplot(contrib[rowSums(abs(contrib)) > 1.5,], is.cor = FALSE)


# Sankey diagram ####

## pathways data ####

# Create data which can be used for Sankey
set.seed(111)

d = genome_paths %>% 
  filter(Completeness == 100) %>% 
  # filter(Name %in% diff_pathways) %>%  ## filter the pathways that seem to be different
  left_join(metadata %>% 
              select(Genome, phylogroup)) %>%
  drop_na(phylogroup) %>% 
  mutate(phylogroup = case_when(phylogroup == 'E or cladeI' ~ 'E',
                                TRUE ~ phylogroup)) %>% 
  select(Name, phylogroup) %>% 
  # filter(Name == 'curcumin degradation')
  filter(Name %in% test_paths)


# Step 1

df <- d %>%
  make_long(Name, phylogroup)
df

# general chart
df %>% 
ggplot(aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node, 
               fill = factor(node), 
               label = node)) + 
  geom_sankey(flow.alpha = 0.5, 
              node.color = "black",
              show.legend = FALSE) +
  geom_sankey_label(size = 3, 
                    color = "black", 
                    fill= "white", hjust = -0.3) +
  theme_bw() +  
  theme(legend.position = "none") + 
  theme(axis.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank())+
  scale_fill_viridis_d(option = "inferno")+
  labs(fill = 'Nodes')

quartz.save(file = '../exploration/chisq_plots/sankey_test.pdf',
            type = 'pdf', height = 6, width = 7)

# highlight nodes
ggplot(df, aes(x = x
               , next_x = next_x
               , node = node
               , next_node = next_node
               , fill = factor(node)
               , label = node)) + 
  geom_sankey(flow.alpha = 0.5, 
              node.color = "black",
              show.legend = FALSE) +
  geom_sankey_label(size = 3, 
                    color = "black", 
                    fill= "white", 
                    hjust = -0.3) +
  theme_bw()+  
  theme(legend.position = "none") + 
  theme(axis.title = element_blank()
        , axis.text.y = element_blank()
        , axis.ticks = element_blank()  
        , panel.grid = element_blank())+
  # scale_fill_viridis_d(option = "inferno")+
  scale_fill_manual(values = c('curcumin degradation' = "red"
                               # , 'B2'  = "#31C9F7"
                               
  ) ) +
  labs(fill = 'Nodes')


quartz.save(file = '../exploration/chisq_plots/sankey_test_highlight.pdf',
            type = 'pdf', height = 6, width = 7)



# create a function that plots a sankey diagram of a selected pathway
plotSankey = function(db, pathway, save.plot = TRUE){

  # filter the database for the pathway
  db_pathway = db %>% 
    filter(Name == pathway) %>% 
    make_long(Name, phylogroup)
  # plot the sankey diagram
  p = ggplot(db_pathway, aes(x = x
                         , next_x = next_x
                         , node = node
                         , next_node = next_node
                         , fill = factor(node)
                         , label = node)) + 
    geom_sankey(flow.alpha = 0.5, 
                node.color = "black",
                show.legend = FALSE) +
    geom_sankey_label(size = 3, 
                      color = "black", 
                      fill= "white", hjust = -0.3) +
    theme_bw() +  
    theme(legend.position = "none") + 
    theme(axis.title = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          panel.grid = element_blank())+
    scale_fill_viridis_d(option = "inferno")+
    labs(fill = 'Nodes')
  if (save.plot == TRUE){
    # save the plot with the glue function
    ggsave(file = glue('../exploration/chisq_plots/sankey_individuals/sankey_{pathway}.pdf'),
           height = 6, width = 7)
  } 
  return(p)
}

db_sankey = genome_paths %>% 
  filter(Completeness == 100) %>% 
  # filter(Name %in% diff_pathways) %>%  ## filter the pathways that seem to be different
  left_join(metadata %>% 
              select(Genome, phylogroup)) %>%
  drop_na(phylogroup) %>% 
  mutate(phylogroup = case_when(phylogroup == 'E or cladeI' ~ 'E',
                                TRUE ~ phylogroup)) 

plotSankey(db_sankey, 'D-erythronate degradation I')

# remove a problematic name
diff_paths_sankey = diff_pathways[diff_pathways != "coenzyme B/coenzyme M regeneration I (methanophenazine-dependent)"]


for (path in diff_paths_sankey){
  cat(glue('Printing {path} pathway\n\n'))
  plotSankey(db_sankey, path)
}


plotSankey(db_sankey, 'L-arabinose degradation II', save.plot = F)

db_pathway = db_sankey %>% 
  filter(Name == 'curcumin degradation') %>% 
  make_long(Name, phylogroup)



# get pathway PA list in a file ####

genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0) %>% 
  left_join(metadata %>% 
              select(Genome, phylogroup, Well, PG)) %>% 
  select(Genome, phylogroup:PG, everything()) %>% 
  pivot_longer(cols = `folate transformations III (E. coli)`:`L-arabinose degradation II`,
               names_to = 'Pathway', values_to = 'Presence') %>% 
  mutate(`Is different?` = case_when(Pathway %in% diff_pathways ~ 'Yes',
                                     TRUE ~ 'No')) %>% 
  write.xlsx('../exploration/tables/pathway_PA_metadata.xlsx')







# glyphosate degradation ####


## assoc Chi2 ####

test_paths = c('ATP biosynthesis', 'isoniazid activation', # present in all
               'Agmatine Transport', 'mannitol cycle', # present in all
               'homotaurine degradation', 'acetoin degradation', # present in all
               'curcumin degradation', 'maltose degradation', # present in a few
               'salmochelin degradation', 'glyphosate degradation III') # present in a few

test_dt = dt[test_paths,]

assoc(test_dt, shade = TRUE, las=1,
      labeling_args = list(
        offset_labels = c(left = -4),
        offset_varnames = c(left = 0),
        rot_labels = c(left = 0)))

quartz.save(file = '../exploration/chisq_plots/glyphosate_deg_assoc_test.pdf',
            type = 'pdf', height = 10, width = 7)



## PA of path in cohort ####

genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0)



length(unique(genome_paths$Genome))

genome_paths %>% 
  filter(str_detect(Name, 'glyphosate')) %>% 
  filter(Completeness > 0) %>% 
  mutate(complete = case_when(Completeness == 100 ~ 'Complete',
                              Completeness < 100 ~ 'Partial'),
         .before = VagueReactions) %>% 
  left_join(metadata, by = c('Genome')) %>% 
  count(complete) %>% 
  # drop_na(Broadphenotype) %>% 
  ggplot(aes(x = complete, y = n, fill =  complete)) +
  geom_histogram(stat= 'identity', show.legend = F, color = 'black') +
  scale_fill_manual(values = c('#339EF5', '#FFDC29')) +
  labs(x = 'Pathway completeness',
       y = 'Count',
       title = 'Glyphosate degradation III') +
  theme_cowplot(15) +
  theme(plot.title = element_text(hjust = 0.5))


quartz.save(file = '../exploration/glyphosate_deg_completeness.pdf',
            type = 'pdf', height = 6, width = 4)




# by bacterial phenotype
genome_paths %>% 
  filter(str_detect(Name, 'glyphosate')) %>% 
  filter(Completeness > 0) %>% 
  mutate(complete = case_when(Completeness == 100 ~ 'Complete',
                              Completeness < 100 ~ 'Partial'),
         .before = VagueReactions) %>% 
  left_join(metadata, by = c('Genome')) %>% 
  drop_na(Broadphenotype) %>% 
  filter(Broadphenotype != 'Unknown') %>% 
  ggplot(aes(complete, fill = Broadphenotype)) +
  geom_bar(stat = 'count', position = 'dodge') +
  facet_wrap(~Broadphenotype, scales = 'free_y') +
  labs(x = 'Pathway completeness',
       y = 'Count',
       title = 'Glyphosate degradation III') +
  theme_cowplot(15) +
  theme(plot.title = element_text(hjust = 0.5))

quartz.save(file = '../exploration/glyphosate_deg_completeness_phenotype.pdf',
            type = 'pdf', height = 6, width = 9)

genome_paths %>% 
  filter(str_detect(Name, 'glyphosate')) %>% 
  filter(Completeness > 0) %>% 
  mutate(complete = case_when(Completeness == 100 ~ 'Complete',
                              Completeness < 100 ~ 'Partial'),
         .before = VagueReactions) %>% 
  left_join(metadata, by = c('Genome')) %>% 
  drop_na(Broadphenotype) %>% 
  filter(Broadphenotype == 'Laboratory strain') %>% 
  select(Name, Prediction, Completeness, Genome, Strainname, Broadphenotype, phylogroup) %>% 
  write_csv('../exploration/glyphosate_lab_strains.csv')





# annotation by Prokka ----------------------------------------------------

# this file hasb been created in the script and R Session from the folder
# /Users/danmarti/Documents/MRC_postdoc/My_projects/pangenome_functions


PG_annotation_prokka = read_csv("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/gapseq/PG_annotation_prokka.csv")


unique(PG_annotation_prokka$genome)

gff_count = PG_annotation_prokka %>% 
  mutate(hypo = case_when(product == 'hypothetical protein' ~ 'hypothetical protein',
                          TRUE ~ 'annotated')) %>% 
  group_by(genome) %>% 
  count(hypo)


gff_count %>% 
  filter(!(genome %in% c('6', 'SPC_3.1', '2.3','OP50'))) %>% 
  # filter(hypo == 'annotated') %>% 
  ggplot(aes(x = n, fill= hypo)) + 
  geom_histogram(bins = 100, binwidth = 50) +
  labs(
    x = 'Genes',
    y = 'Genome count'
  ) +
  scale_fill_manual(values = c('#FA9228', '#8558F0')) +
  theme_cowplot(15) +
  theme(legend.position = 'top',
        legend.justification = "center") +
  guides(fill = guide_legend(title='Category'))

ggsave(here('../exploration', 'PG_annotation', 'annotated_genes_histogram.pdf'),
       height = 7, width = 9)


gff_count %>% 
  filter(!(genome %in% c('6', 'SPC_3.1', '2.3','OP50'))) %>% 
  ggplot(aes(y = fct_reorder(genome, n, .desc = T), x = n, fill = hypo)) +
  geom_bar(position = 'stack', stat = 'identity', width = 1) +
  labs(
    x = 'Number of genes',
    y = 'Genomes'
  ) +
  scale_fill_manual(values = c('#FA9228', '#8558F0'),
                    labels = c('Annotated proteins',
                               'Hypothetical protein')) +
  scale_y_discrete(limits=rev) +
  # scale_x_discrete(limits=rev) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(here('../exploration', 'PG_annotation', 'annotated_genes_barplot.pdf'),
       height = 7, width = 9)


# most repeated genes

PG_annotation_prokka %>% 
  mutate(genome = str_sub(genome, end = -5)) %>% 
  drop_na(gene) %>% 
  filter(product != 'hypothetical protein') %>% 
  count(gene, sort = TRUE) %>% 
  ggplot(aes(x = n)) +
  geom_bar(position = 'identity', width = 1) +
  labs(
    x = 'Genome count',
    y = 'Gene count'
  ) +
  theme_cowplot(14)

ggsave(here('../exploration', 'PG_annotation', 'annotated_genes_barplot.pdf'),
       height = 7, width = 9)

# gff_db %>% 
#   dplyr::select(-value) %>% 
#   mutate(genome = str_sub(genome, end = -5)) %>% 
#   drop_na(gene) %>% 
#   filter(product != 'hypothetical protein') %>% 
#   separate(gene, into = c('gene', 'trash'), sep = '_') %>% 
#   count(gene, sort = TRUE) %>% 
#   ggplot(aes(x = n)) +
#   geom_bar(position = 'identity', width = 1) +
#   labs(
#     x = 'Genome count',
#     y = 'Gene count'
#   ) +
#   theme_cowplot(14)
# 
# 
# 
# cosa = gff_db %>% 
#   dplyr::select(-value) %>% 
#   mutate(genome = str_sub(genome, end = -5)) %>% 
#   drop_na(gene) %>% 
#   filter(product != 'hypothetical protein') %>% 
#   separate(gene, into = c('gene', 'trash'), sep = '_') %>% 
#   count(gene, sort = TRUE) 







