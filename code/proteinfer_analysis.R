
# libraries ---------------------------------------------------------------


library(tidyverse)
library(readr)
library(glue)
library(cowplot)
library(here)
library(ggrepel)

theme_set(theme_cowplot(14))


two_cols = c('#EB7D0E', '#0EB7EB')
three_cols = c('#A428EB', '#EB9D11', '#28EB83')
four_cols = c('#EB9C1A', '#1081EB', '#EB21C4', '#34C90A')
five_cols = c('#02C9F7', '#24D402', '#EBB70D', '#D42402', '#8C07F5')

dir.create('exploration/proteinfer', showWarnings = F)

# read data ---------------------------------------------------------------

proteinfer = read_delim("raw_data/proteinfer/combined_file_fold2.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) %>% 
  rename(gene = sequence_name) %>% 
  drop_na(gene)

protein_lengths = read_csv("raw_data/proteinfer/protein_lengths.csv") %>% 
  select(-`...1`) %>% 
  rename(gene = name)

proteinfer_functions = proteinfer %>% 
  distinct(predicted_label, description)

# exploration -------------------------------------------------------------

length(unique(proteinfer$gene)) / length(gene_PA$Gene)

# predicted genes vs not predicted
# build a df with every gene

genes_not_proteinfer = setdiff(gene_PA$Gene, unique(proteinfer$gene))

proteinfer_full = tibble(gene = genes_not_proteinfer) %>% 
  full_join(proteinfer) %>% 
  mutate(annotated = case_when(is.na(predicted_label) ~ 'No',
                               TRUE ~ 'Yes'))

proteinfer_full %>% 
  write_csv('tables/proteinfer_full.csv')


proteinfer_full %>% 
  distinct(gene, .keep_all = T) %>%
  count(annotated) %>% 
  mutate(prop = (n / sum(n)) * 100) %>% 
  ggplot(aes(x = annotated, y = prop, fill = annotated)) + 
  labs(x = 'Predicted by ProteInfer?',
       y = "%") +
  scale_fill_manual(values = two_cols) +
  geom_col(color = 'black', show.legend = F)

ggsave(here('exploration', 'proteinfer', 
            'predicted_prop.pdf'), height = 6, width = 5)

proteinfer_full %>% 
  distinct(gene, .keep_all = T) %>%
  mutate(group = ifelse(str_detect(gene, 'group'), 'Group', 'Gene name')) %>% 
  count(annotated, group) %>% 
  drop_na() %>% 
  pivot_wider(names_from = group, values_from = n) %>% 
  select(-annotated) %>% 
  chisq.test()

proteinfer_full %>% 
  distinct(gene, .keep_all = T) %>%
  mutate(group = ifelse(str_detect(gene, 'group'), 'Group', 'Gene name')) %>% 
  count(annotated, group) %>% 
  drop_na() %>% 
  mutate(prop = (n / sum(n)) * 100) %>% 
  ggplot(aes(x = annotated, y = prop, fill = group)) + 
  labs(x = 'Predicted by ProteInfer?',
       y = "%") +
  scale_fill_manual(values = two_cols) +
  geom_col(color = 'black') +
  annotate('text', x = 1, y = 75, 
           label = 'Chisq test\n p-value < 2.2e-16',
           size = 5) +
  theme(legend.position = c(0.2, 0.6))

ggsave(here('exploration', 'proteinfer', 
            'predicted_group_prop.pdf'), height = 6, width = 5)

# TODO: find a better representation for unnanotated proteins by prokka or 
# by COGs

# predictions confidence density plot
proteinfer_full %>% 
  filter(annotated == 'Yes') %>% 
  ggplot(aes(confidence)) +
  geom_density( fill="lightblue")

ggsave(here('exploration', 'proteinfer', 
            'confidence_density.pdf'), height = 4, width = 5)


## GO terms ----- 
# most annotated GO terms
go_prot = proteinfer %>% 
  count(predicted_label) %>% 
  drop_na() %>% 
  filter(str_detect(predicted_label, 'GO:')) %>% 
  arrange(desc(n)) %>% 
  left_join(proteinfer_functions)

go_prot%>% 
  head(500) %>% 
  mutate(description = ifelse(n > 2886, description, '')) %>% 
  ggplot(aes(x = fct_reorder(predicted_label, n, .desc = T), y = n)) +
  geom_col(fill = 'black') +
  geom_text_repel(aes(label = description), box.padding = 0.5,
                  max.overlaps = Inf) +
  labs(x = 'GO terms (predicted)') +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(here('exploration', 'proteinfer', 
            'GO_terms_histogram.pdf'), height = 6, width = 8)


## Pfam -----

pfam_prot = proteinfer %>% 
  count(predicted_label) %>% 
  drop_na() %>% 
  filter(str_detect(predicted_label, 'Pfam:')) %>% 
  arrange(desc(n)) %>% 
  left_join(proteinfer_functions)

pfam_prot %>% 
  head(500) %>% 
  mutate(description = ifelse(n > 70, description, '')) %>% 
  ggplot(aes(x = fct_reorder(predicted_label, n, .desc = T), y = n)) +
  geom_col(fill = 'black') +
  geom_text_repel(aes(label = description), box.padding = 0.5,
                  max.overlaps = Inf) +
  labs(x = 'Pfam terms (predicted)') +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(here('exploration', 'proteinfer', 
            'Pfam_terms_histogram.pdf'), height = 6, width = 8)


## EC numbers -----

ec_prot = proteinfer %>% 
  count(predicted_label) %>% 
  drop_na() %>% 
  filter(str_detect(predicted_label, 'EC:')) %>% 
  arrange(desc(n)) %>% 
  left_join(proteinfer_functions)

ec_prot %>% 
  head(500) %>% 
  unite(labels, predicted_label, description, remove=F, sep = ": ") %>% 
  mutate(labels = ifelse(n > 136, labels, '')) %>% 
  ggplot(aes(x = fct_reorder(predicted_label, n, .desc = T), y = n)) +
  geom_col(fill = 'black') +
  geom_text_repel(aes(label = labels), box.padding = 0.6,
                  max.overlaps = Inf) +
  labs(x = 'EC terms (predicted)') +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(here('exploration', 'proteinfer', 
            'EC_terms_histogram.pdf'), height = 8, width = 12)





proteinfer_full %>% 
  drop_na(predicted_label) %>%
  select(gene, predicted_label) %>% 
  count(gene) %>%
  ggplot(aes(x = fct_reorder(gene, n), y = n)) + 
  geom_bar(stat = 'identity', fill = 'black') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


proteinfer_full %>% 
  drop_na(predicted_label) %>%
  select(gene, predicted_label) %>% 
  count(gene) %>%
  ggplot(aes(n)) +
  geom_histogram(bins = 70, fill = 'black') +
  theme()


# correlation  ------------------------------------------------------------



proteinfer_full %>% 
  drop_na(predicted_label) %>% 
  group_by(gene) %>% 
  count() %>% 
  left_join(protein_lengths) %>% 
  ggplot(aes(x = n, y = length)) +
  geom_point(alpha = 0.1, shape = 21, size = 2 , fill ='darkred') +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 6000,
                   size = 6) +
  labs(
    x = 'Number of categories predicted',
    y = 'Protein length (in aa)'
  )


ggsave(here('exploration', 'proteinfer', 
            'correlation_cats_length.pdf'), height = 6, width = 8)



proteinfer_full %>% 
  drop_na(predicted_label) %>% 
  group_by(gene) %>% 
  count() %>% 
  left_join(protein_lengths) %>% 
  mutate(groups = case_when(str_detect(gene, 'group') ~ 'group',
                            TRUE ~ 'gene')) %>% 
  ggplot(aes(x = n, y = length, color = groups)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = "pearson", label.x = 50, 
                   size = 6) +
  scale_color_manual(values = two_cols) +
  labs(
    x = 'Number of categories predicted',
    y = 'Protein length (in aa)'
  )

ggsave(here('exploration', 'proteinfer', 
            'correlation_cats_length_groups.pdf'), height = 6, width = 8)


# Multivariate ----------------------------------------------------------------

proteinfer_pca_ready = proteinfer_full %>% 
  drop_na(predicted_label) %>%
  select(gene, predicted_label) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = predicted_label, 
              values_from = presence, 
              values_fill = 0)

proteinfer_pca_ready %>% 
  write_csv(here('exploration', 'proteinfer',
                 'proteinfer_categories_wide.csv'))



# size of a K-12 strain: MG1655 ----------------------------

k12_gene_sizes = gene_PA_long %>% 
  filter(Genome == 'NT12001_189', presence == 1) %>% 
  left_join(protein_lengths %>% 
              rename(Gene = gene)) %>% 
  arrange(desc(length))


gene_PA_long %>% 
  filter(Genome == 'NT12009_154', presence == 1) %>% 
  left_join(protein_lengths %>% 
              rename(Gene = gene)) %>% 
  arrange(desc(length))

gene_max_length = gene_PA_long %>% 
  filter(presence == 1) %>% 
  left_join(protein_lengths %>% 
              rename(Gene = gene)) %>% 
  arrange(desc(length)) %>% 
  group_by(Genome) %>% 
  slice_head(n = 1)

# dotplot of gene max length by phylogroup
gene_max_length %>% 
  left_join(metadata %>% 
              select(Genome, phylogroup, Broadphenotype)) %>% 
  ggplot(aes(x = fct_reorder(Genome, length),
             y = length,
             color = phylogroup)) +
  geom_point(size = 4) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(here('exploration', 'proteinfer',
            'max_gene_length_dots.pdf'),
       height = 5, width = 7)

# boxplot of gene max length by phylogroup
gene_max_length %>% 
  left_join(metadata %>% 
              select(Genome, phylogroup, Broadphenotype)) %>% 
  filter(phylogroup != 'E or cladeI') %>% 
  drop_na(phylogroup) %>% 
  ggplot(aes(x = phylogroup, y = length, fill = phylogroup)) +
  geom_boxplot(show.legend = F) +
  geom_point(position = position_jitterdodge(), 
             alpha = 0.2,
             show.legend = F) +
  labs(y = 'Length of largest gene in genome',
       x = NULL)

ggsave(here('exploration', 'proteinfer',
            'max_gene_length_boxplot.pdf'),
       height = 5, width = 7)


# gene length density

gene_length.sum = gene_PA_long %>% 
  filter(presence == 1) %>% 
  left_join(metadata) %>% 
  drop_na(phylogroup) %>% 
  filter(phylogroup != 'E or cladeI') %>% 
  left_join(protein_lengths %>% 
              rename(Gene = gene)) %>% 
  group_by(phylogroup, Genome) %>% 
  summarise(gene_mean = mean(length, na.rm = TRUE),
            gene_median = median(length, na.rm = TRUE))


gene_length.sum %>% 
  ggplot(aes(x = phylogroup, y = gene_mean, fill = phylogroup)) +
  geom_boxplot(show.legend = F) +
  geom_point(position = position_jitterdodge(), 
             alpha = 0.2,
             show.legend = F) 
  

# interpro results --------------------------

interpro = read_delim("raw_data/interpro/interproscan_PG.tsv", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE,
                      trim_ws = TRUE)

interpro = interpro %>% 
  rename(gene = X1,
         seq_length = X3,
         analysis = X4, 
         signature_accession = X5,
         signature_description = X6, 
         start = X7,
         end = X8,
         score = X9,
         status = X10,
         date =  X11, 
         interpro_accession = X12,
         interpro_description = X13,
         go_terms = X14) %>% 
  select(-X2)

# how many genes have a true GO prediction
interpro %>% 
  drop_na(go_terms) %>% 
  filter(go_terms != '-') %>% 
  distinct(gene) %>% 
  count


# filter the genes that have GO prediction and separate rows, keep only distinct
# GO terms for each gene as in Proteinfer
interpro_go = interpro %>% 
  drop_na(go_terms) %>% 
  filter(go_terms != '-') %>% 
  separate_rows(go_terms, sep = '\\|') %>% 
  distinct(gene, go_terms, .keep_all = TRUE)


# are the genes from interpro in the proteinfer list?
length(setdiff(interpro_go$gene, proteinfer$gene))



# let's make a barplot of the missing 
interpro_go_genes = unique(interpro_go$gene)

gene_PA %>% 
  select(Gene) %>% 
  mutate(interpro_predict = case_when(Gene %in% interpro_go_genes ~ 'Yes',
                                      TRUE ~ 'No'),
         annotated = case_when(str_detect(Gene, 'group') ~ 'group',
                               TRUE ~ 'annotated')) %>%
  group_by(interpro_predict, annotated) %>% 
  mutate(interpro_predict = factor(interpro_predict),
         annotated = factor(annotated)) %>% 
  count %>% 
  pivot_wider(values_from = n, names_from = annotated) %>% 
  ungroup() %>% 
  select(-interpro_predict) %>% 
  chisq.test()
  

gene_PA %>% 
  select(Gene) %>% 
  mutate(interpro_predict = case_when(Gene %in% interpro_go_genes ~ 'Yes',
                                      TRUE ~ 'No'),
         group = case_when(str_detect(Gene, 'group') ~ 'Group',
                           TRUE ~ 'Gene name')) %>%
  group_by(interpro_predict, group) %>% 
  count %>% 
  ungroup %>% 
  mutate(prop = n/sum(n) * 100) %>% 
  ggplot(aes(x = interpro_predict, y = prop, fill = group)) + 
  labs(x = 'Predicted by InterPro?',
       y = "%") +
  scale_fill_manual(values = two_cols) +
  geom_col(color = 'black') +
  theme(legend.position = c(0.6, 0.67)) +
  ylim(0, 80) +
  annotate('text', x = 2, y = 70, label = 'Chisq test\n p-value < 2.2e-16')

ggsave(here('exploration', 'proteinfer', 
            'predicted_group_prop_interpro.pdf'), height = 6, width = 5)



## correlation of GO terms between interpro and proteinfer


corr_go = interpro_go %>% 
  filter(gene %in% intersect(interpro_go$gene, proteinfer$gene)) %>%
  select(gene, go_terms) %>% 
  mutate(method = 'interpro') %>% 
  bind_rows(proteinfer %>% 
              filter(gene %in% intersect(interpro_go$gene, proteinfer$gene),
                     str_detect(predicted_label, 'GO:')) %>% 
              select(gene, go_terms = predicted_label) %>% 
              mutate(method = 'proteinfer')) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = method, values_from = presence) %>% 
  mutate_all(~replace_na(.,0))

corr_go %>% 
  ggplot(aes(x = interpro, y = proteinfer)) +
  geom_point() +
  geom_smooth(method = 'lm')
  

cor(corr_go$interpro, corr_go$proteinfer)

# lollipop plot for the elements in each method
corr_go %>% 
  mutate(presence = case_when(proteinfer == 1 & interpro == 1 ~ 'Both',
                              proteinfer == 1 & interpro == 0 ~ 'ProteInfer',
                              proteinfer == 0 & interpro == 1 ~ 'InterPro')) %>% 
  count(presence) %>% 
  mutate(log_n = log2(n)) %>% 
  ggplot(aes(x = log_n, y = presence)) +
  geom_segment(aes(y=presence, yend=presence, x=0, xend=log_n)) +
  geom_point(size=5, color="red", fill=alpha("orange", 0.3), 
             alpha=0.7, shape=21, stroke=2) +
  # xlim(0, 230000) +
  labs(
    x = 'Number of elements (in log2)',
    y = NULL
  ) +
  ggrepel::geom_text_repel(aes(label = round(n)),
                           box.padding = 0.2,
                           nudge_x = .25, 
                           nudge_y = 0.19, 
                           size = 5 )
ggsave(here('exploration', 'proteinfer', 
            'interpro_vs_proteinfer_lollipop.pdf'), height = 5, width = 6.5)


corr_go %>% 
  filter(proteinfer == 1 & interpro == 0)
corr_go %>% 
  filter(proteinfer == 0 & interpro == 1)
corr_go %>% 
  filter(proteinfer == 1 & interpro == 1)


# number of categories per method and per gene
interpro_go %>% 
  filter(gene %in% intersect(interpro_go$gene, proteinfer$gene)) %>%
  select(gene, go_terms) %>% 
  mutate(method = 'interpro') %>% 
  bind_rows(proteinfer %>% 
              filter(gene %in% intersect(interpro_go$gene, proteinfer$gene),
                     str_detect(predicted_label, 'GO:')) %>% 
              select(gene, go_terms = predicted_label) %>% 
              mutate(method = 'proteinfer')) %>% 
  mutate(presence = 1) %>% 
  group_by(gene, method) %>% 
  count() %>% 
  group_by(method) %>% 
  summarise(Mean = mean(n, na.rm = T),
            Median = median(n, na.rm = T),
            SD = sd(n, na.rm = T))



# drug related proteins ---------------------------------------------------


go_prot %>% 
  filter(str_detect(description, 'drug'))


proteinfer %>% 
  filter(str_detect(description, 'drug metabolic process')) %>% view

library(FGNet)
plotGoAncestors(c("GO:0051603", "GO:0019941", "GO:0051128","GO:0044265"),
                plotOutput="dynamic")
plotGoAncestors(c("GO:0000152","GO:0043234", "GO:0044446", "GO:0043227"))








