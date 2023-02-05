
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

proteinfer %>% 
  count(predicted_label) %>% 
  drop_na() %>% 
  filter(str_detect(predicted_label, 'GO:')) %>% 
  arrange(desc(n)) %>% 
  head(1000) %>% 
  left_join(proteinfer_functions) %>% view

proteinfer %>% 
  count(predicted_label) %>% 
  drop_na() %>% 
  filter(str_detect(predicted_label, 'GO:')) %>% 
  arrange(desc(n)) %>% 
  head(500) %>% 
  left_join(proteinfer_functions) %>% 
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

proteinfer %>% 
  count(predicted_label) %>% 
  drop_na() %>% 
  filter(str_detect(predicted_label, 'Pfam:')) %>% 
  arrange(desc(n)) %>% 
  head(500) %>% 
  left_join(proteinfer_functions) %>% 
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


## Pfam -----

proteinfer %>% 
  count(predicted_label) %>% 
  drop_na() %>% 
  filter(str_detect(predicted_label, 'EC:')) %>% 
  arrange(desc(n)) %>% 
  head(500) %>% 
  left_join(proteinfer_functions) %>% 
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




