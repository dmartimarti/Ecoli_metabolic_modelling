
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

proteinfer = read_delim("raw_data/proteinfer/combined_file.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) %>% 
  rename(gene = sequence_name)

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
  mutate(prop = (n / sum(n)) * 100) %>% 
  ggplot(aes(x = annotated, y = prop, fill = group)) + 
  labs(x = 'Predicted by ProteInfer?',
       y = "%") +
  scale_fill_manual(values = two_cols) +
  geom_col(color = 'black') +
  theme(legend.position = c(0.2, 0.6))

ggsave(here('exploration', 'proteinfer', 
            'predicted_group_prop.pdf'), height = 6, width = 5)

# TODO: find a better representation for unnanotated proteins by prokka or 
# by COGs


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



# drug related proteins ---------------------------------------------------


go_prot %>% 
  filter(str_detect(description, 'drug'))


proteinfer %>% 
  filter(str_detect(description, 'drug metabolic process')) %>% view

library(FGNet)
plotGoAncestors(c("GO:0051603", "GO:0019941", "GO:0051128","GO:0044265"),
                plotOutput="dynamic")
plotGoAncestors(c("GO:0000152","GO:0043234", "GO:0044446", "GO:0043227"))


