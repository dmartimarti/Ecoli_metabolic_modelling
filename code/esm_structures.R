library(tidyverse)
library(NGLVieweR)
library(glue)
library(bio3d)
library(cowplot)
library(ggpubr)

theme_set(theme_cowplot(14))



# pdb viewer --------------------------------------------------------------


root = '/Users/danmarti/Documents/MRC_postdoc/My_projects/Ecoli_metabolic_modelling/prot_structures_esm/pdb_files'
proteins = list.files(root, pattern = 'pdb')

prot_df = proteins %>% 
  as_tibble() %>% 
  mutate(id = 1:n()) %>% 
  rename(protein = value)

list.files('PG_protein_vis/data', pattern = 'pdb') %>% 
  write('PG_protein_vis/data/prot_list.txt')

plot_protein <- function(prot_name = 'apaH') { 
  prot_name = glue("{prot_name}.pdb")
  prot_id = prot_df %>% 
    filter(protein == prot_name) %>% 
    head(1) %>% 
    pull(id)
  
  print(glue("Showing protein {prot_name}"))
  NGLVieweR(glue('{root}/{proteins[prot_id]}')) %>%
    stageParameters(backgroundColor = "white") %>% 
    addRepresentation("cartoon",
      param = list(name = "cartoon", colorScheme = "bfactor") 
    )
}


plot_protein('yiaN_1~~~yiaN_2~~~yiaN_3')




# reading pdb outputs -----------------------------------------------------


pbs_metadata = tibble()
for (out in list.files('pbs_outputs', pattern = 'esmfold_batch')) {
  temp = read_delim(glue("pbs_outputs/{out}"), 
             delim = "|", escape_double = FALSE, col_names = FALSE, 
             trim_ws = TRUE, skip = 6) %>% 
    filter(!str_detect(`X4`, 'CUDA out of memory'))
  
  pbs_metadata = pbs_metadata %>% 
    bind_rows(temp)
}

pbs_metadata %>%
  write_csv('exploration/esm_metadata.csv')


# model performance -------------------------------------------------------



pbs_metadata = pbs_metadata %>% 
  mutate(protein = str_extract(`X4`, "for\\s(.*?)\\s"),
         length = str_extract(`X4`, "length\\s(\\d+)"),
         pLDDT = str_extract(`X4`, "pLDDT\\s(\\d+\\.\\d+)"),
         pTM = str_extract(`X4`, "pTM\\s(\\d+\\.\\d+)")) %>% 
  mutate(protein = str_sub(protein, start = 5, end = -2),
         length = str_sub(length, start = 7),
         pLDDT = str_sub(pLDDT, start =6), 
         pTM = str_sub(pTM, start = 4)) %>% 
  mutate_at(vars(length, pLDDT, pTM), as.numeric) %>% 
  select(protein:pTM) %>% 
  mutate(plddt_class = case_when(pLDDT >= 90 ~ 'Very high',
                                 pLDDT < 90 & pLDDT >= 70 ~ 'Confident',
                                 pLDDT < 70 & pLDDT >= 50 ~ 'Low',
                                 pLDDT < 50 ~ 'Very low'),
         plddt_class = factor(plddt_class, 
                              levels = c('Very high', 'Confident',
                                        'Low', 'Very low'))) %>% 
  left_join(prot_df)


plddt_cols = c('#21A9E0','#FAEA2E','#2F40D9','#E08F42')
# correlation between length and structure quality
pbs_metadata %>% 
  ggplot(aes(x = pLDDT, y = length, color = plddt_class)) +
  geom_point(alpha = .3) +
  geom_smooth(method = 'lm', color='black') +
  scale_color_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  stat_cor(method = "pearson", label.x = 30, label.y = 1250, color = 'black',
           size = 6) +
  labs(y = 'Protein length (aa)',
       color = 'Model mean\nconfidence')
ggsave('exploration/length_plddt.pdf', height = 7, width = 9)


mean_dist = round(mean(pbs_metadata$pLDDT),2)
median_dist = round(median(pbs_metadata$pLDDT),2)
pbs_metadata %>% 
  ggplot(aes(x = pLDDT, fill = plddt_class, color = plddt_class)) +
  geom_histogram(binwidth=1, width = 0.1) +
  geom_vline(xintercept = mean_dist, color = 'blue') +
  annotate("text", x = mean_dist - 5, 
           y = 850, label = glue('mean: \n{mean_dist}'), color = 'blue') +
  geom_vline(xintercept = median_dist, color = 'red') +
  annotate("text", x = median_dist + 5, 
           y = 850, label = glue('median: \n{median_dist}'), color = 'red') +
  scale_fill_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  scale_color_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  labs(fill = 'Model mean\nconfidence',
       color = 'Model mean\nconfidence',
       y = 'count')

ggsave('exploration/plddt_histo.pdf', height = 7, width = 9)

pbs_metadata %>% 
  group_by(plddt_class) %>% 
  count() %>% 
  ggplot(aes(x = plddt_class, y=n , fill = plddt_class, color = plddt_class)) +
  geom_col(show.legend = F) +
  scale_fill_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  scale_color_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  labs(fill = 'Model mean\nconfidence',
       color = 'Model mean\nconfidence',
       y = 'count', x = 'Model mean confidence')
ggsave('exploration/plddt_cases.pdf', height = 5, width = 7)



mean_dist = round(mean(pbs_metadata$length),2)
median_dist = round(median(pbs_metadata$length),2)
pbs_metadata %>% 
  ggplot(aes(x = length, fill= plddt_class, color = plddt_class)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = mean_dist, color = 'blue') +
  annotate("text", x = mean_dist + 55,
           y = 1280, label = glue('mean: \n{mean_dist}'), color = 'blue') +
  geom_vline(xintercept = median_dist, color = 'red') +
  annotate("text", x = median_dist - 55,
  y = 1280, label = glue('median: \n{median_dist}'), color = 'red') +
  scale_fill_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  scale_color_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  labs(fill = 'Model mean\nconfidence',
       color = 'Model mean\nconfidence',
       y = 'count') +
  theme(legend.position = c(0.5, 0.5))
ggsave('exploration/prot_length_histo.pdf', height = 7, width = 8)


medians = pbs_metadata %>% 
  group_by(plddt_class) %>% 
  summarise(median_length = median(length))

pbs_metadata %>% 
  ggplot(aes(x = length, fill= plddt_class, color = plddt_class)) +
  geom_histogram(bins = 100, show.legend = F, width=1) +
  scale_fill_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  scale_color_manual(values = c('#2F40D9','#21A9E0','#FAEA2E','#E08F42')) +
  labs(fill = 'Model mean\nconfidence',
       y = 'count') +
  geom_vline(data = medians, aes(xintercept = median_length)) +
  facet_wrap(~plddt_class, scales = 'free_y', ncol = 1) 

ggsave('exploration/prot_length_histo_facet.pdf', height = 9, width = 4)


# stats -------------------------------------------------------------------

## pLDDT vs annotated genes -----

# annotated genes have better quality? 
group_df = pbs_metadata %>% 
  mutate(known = case_when(str_detect(protein, 'group') ~ 'group',
                           TRUE ~ 'gene'))
# there's statistical diff between the groups and annotated genes
model = lm(pLDDT ~ known, data = group_df)
summary(model)

summary(model)$coeff[2,1]

pbs_metadata %>% 
  mutate(known = case_when(str_detect(protein, 'group') ~ 'group',
                           TRUE ~ 'gene')) %>% 
  ggplot(aes(x = known, y = pLDDT, fill = known)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#21A9E0','#E08F42')) +
  annotate('text', x = 1.5, y = 100, label = 'p-val<2e-16', size = 5) +
  annotate('text', x = 1.5, y = 95, label = glue('estimate: {round(summary(model)$coeff[2,1],2)}'), 
           size = 4) +
  geom_point(position = position_jitterdodge(), alpha = .1)

ggsave('exploration/plddt_geneVSgroup.pdf', width = 8, height = 6)


## length vs annotated genes -----


# there's statistical diff between the groups and annotated genes
model = lm(length ~ known, data = group_df)
summary(model)

summary(model)$coeff[2,1]

pbs_metadata %>% 
  mutate(known = case_when(str_detect(protein, 'group') ~ 'group',
                           TRUE ~ 'gene')) %>% 
  ggplot(aes(x = known, y = length, fill = known)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#21A9E0','#E08F42')) +
  annotate('text', x = 1.5, y = 1000, label = 'p-val<2e-16', size = 5) +
  annotate('text', x = 1.5, y = 920, label = glue('estimate: {round(summary(model)$coeff[2,1],2)}'), 
           size = 4) +
  geom_point(position = position_jitterdodge(), alpha = .1)

ggsave('exploration/length_geneVSgroup.pdf', width = 8, height = 6)

