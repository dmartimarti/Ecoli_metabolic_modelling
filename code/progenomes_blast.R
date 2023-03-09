library(tidyverse)
library(readr)
library(taxonomizr)
library(cowplot)
library(glue)


theme_set(theme_cowplot(15))

blastx = read_delim("reference_PG_progenomes3.tsv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)


filtered_headers <- read_delim("filtered_headers.tsv", 
                               delim = "\t", escape_double = FALSE, 
                               col_names = FALSE, trim_ws = TRUE) %>% 
  rename(sseqid = X1,
         annotation_project_id = X2,
         tax_id_sample_id = X3,
         sample_id = X4,
         tax_id = X5,
         cog = X6,
         product = X7)


blastx %>% 
  left_join(filtered_headers) %>% 
  filter(evalue < 1e-40) %>% 
  # filter(str_detect(qseqid, 'group')) %>% 
  view

blastx %>% distinct(qseqid) %>% count()

blastx %>% 
  left_join(filtered_headers) %>% 
  drop_na(product) %>% 
  group_by(qseqid) %>% 
  slice_head(n = 4) 

full_blast = blastx %>% 
  left_join(filtered_headers) %>% 
  drop_na(product) %>% 
  group_by(qseqid) 


full_blast %>% 
  write_csv('progenomes_blast_full.csv')

# exploration -------------------------------------------------------------


full_blast %>% 
  group_by(qseqid) %>% 
  arrange(desc(pident)) %>% 
  slice_head(n=1) %>% 
  ungroup %>% 
  ggplot(aes(pident)) +
  annotate(geom = 'rect', xmin = 80, xmax = 101, ymin = 0, ymax = 7600, 
                     fill = 'darkblue', alpha = 0.5) +
  annotate(geom = 'rect', xmin = 20, xmax = 50, ymin = 0, ymax = 7600, 
           fill = 'darkgreen', alpha = 0.5) +
  geom_histogram(bins = 100, color = 'black', fill = 'black') +
  annotate('text', y = 6500, x = 90, label = '~70% of data') +
  annotate('text', y = 6500, x = 35, label = 'Very interesting \nfraction!') +
  labs(
    x = 'Identity to best protein in blast (%)',
    y = 'count'
  )

ggsave('exploration/pident_histogram.pdf', height = 7, width = 8)

# proportion of pident 
full_blast %>% 
  select(qseqid, pident) %>% 
  group_by(qseqid) %>% 
  arrange(desc(pident)) %>% 
  slice_head(n=1) %>% 
  ungroup %>% 
  mutate(group = case_when(pident > 79 ~ 'UP',
                           TRUE ~ 'DOWN')) %>% 
  count(group) %>% 
  mutate(prop =  n / sum(n))


# proportion of pident 
rare_proteins = full_blast %>% 
  select(qseqid, pident) %>% 
  group_by(qseqid) %>% 
  arrange(desc(pident)) %>% 
  slice_head(n=1) %>% 
  filter(pident < 50)




# read data from proteinfer and ESM ---------------------------------------


esm_metadata = read_csv("../tables/esm_metadata.csv")
proteinfer_full = read_csv("../tables/proteinfer_full.csv")
full_blast


full_blast %>% 
  distinct(qseqid)

esm_metadata %>% 
  distinct(protein)

proteinfer_full %>% 
  distinct(gene)

proteinfer_full %>% 
  left_join(esm_metadata %>% 
              rename(gene = protein)) %>% 
  left_join(full_blast %>% 
              ungroup %>% 
              mutate(gene = qseqid) %>% 
              select(-length))


# rare proteins -----------------------------------------------------------

# proportion of pident 
rare_proteins = full_blast %>% 
  select(qseqid, pident) %>% 
  group_by(qseqid) %>% 
  arrange(desc(pident)) %>% 
  slice_head(n=1) %>% 
  filter(pident < 50)

full_blast %>% 
  select(qseqid, pident) %>% 
  group_by(qseqid) %>% 
  arrange(desc(pident)) %>% 
  slice_head(n=1) %>% 
  filter(pident > 50) %>% 
  ungroup %>% 
  count()

# how many well predicted groups do we have? 
rare_prot_func = rare_proteins %>% 
  ungroup %>% 
  # filter(str_detect(qseqid, 'group')) %>% 
  left_join(proteinfer_full %>% 
              mutate(qseqid = gene),
            multiple = "all") %>% 
  left_join(esm_metadata %>% 
              rename(gene = protein))


rare_prot_func %>% count(qseqid)

rare_prot_func %>% 
  count(qseqid) %>% 
  filter(n > 1)

rare_prot_func %>%
  drop_na(predicted_label) %>% 
  group_by(gene, pident) %>% 
  count() %>% 
  arrange(desc(n)) %>% view

# plot histogram of number of terms predicted

rare_prot_func %>% 
  count(qseqid) %>% 
  filter(n > 1) %>% 
  summarise(Med = mean(n))

rare_prot_func %>% 
  count(qseqid) %>%
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1, color = 'black', fill = '#11A9FA') +
  labs(x = 'Categories predicted (proteinfer)',
       caption = 'Rare proteins only (<50% of identity to progenomes DB)')

ggsave('../exploration/proteinfer/rare_prot_histogram.pdf',
       height = 5, width = 6)

rare_prot_func %>% 
  count(qseqid) %>%
  filter(n > 1) %>% 
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1, color = 'black', fill = '#11A9FA') +
  labs(x = 'Categories predicted (proteinfer)',
       caption = 'Rare proteins only (<50% of identity to progenomes DB). count > 1')

ggsave('../exploration/proteinfer/rare_prot_histogram_no1.pdf',
       height = 5, width = 6)

## EC numbers ----------

rare_prot_func %>% 
  distinct(qseqid, pLDDT, .keep_all = TRUE) %>%
  filter(str_detect(predicted_label, 'EC')) %>% 
  mutate(EC_1 = str_sub(predicted_label, start = 4, end = -7),
         .before = confidence)

rare_prot_func %>% 
  distinct(qseqid, pLDDT, .keep_all = TRUE) %>%
  filter(str_detect(predicted_label, 'EC')) %>% 
  count(predicted_label, description) %>% 
  ggplot(aes(n, fct_reorder(predicted_label,n))) +
  geom_col(fill = 'darkgreen') +
  xlim(0, 50) +
  ggrepel::geom_label_repel(aes(label = description),
                            nudge_y = -0.5, nudge_x = 14, 
                            size = 3.5) +
  labs(y = 'EC')

ggsave('../exploration/proteinfer/rare_prot_EC.pdf',
       height = 6, width = 7)

## Pfam domains ----------

rare_prot_func %>% 
  distinct(qseqid, pLDDT, .keep_all = TRUE) %>%
  filter(str_detect(predicted_label, 'Pfam')) %>% 
  mutate(pfam = str_sub(predicted_label, start =6),
         .before = confidence) %>% 
  count(pfam, description) %>% 
  arrange(desc( n)) %>% 
  head(20) %>% 
  ggplot(aes(n, fct_reorder(pfam,n))) +
  geom_col(fill = 'darkred') +
  # xlim(0, 50) +
  ggrepel::geom_label_repel(aes(label = description),
                            nudge_y = -0.5, nudge_x = 14, 
                            size = 3.5) +
  labs(y = 'Pfam',
       caption = 'Top 20 elements')

ggsave('../exploration/proteinfer/rare_prot_Pfam.pdf',
       height = 6, width = 7)


## GO terms ----------

rare_prot_func %>% 
  distinct(qseqid, pLDDT, .keep_all = TRUE) %>%
  filter(str_detect(predicted_label, 'GO')) %>% 
  mutate(go = str_sub(predicted_label, start =1),
         .before = confidence) %>% 
  count(go, description) %>% 
  arrange(desc( n)) %>% 
  filter(n < 120) %>% 
  head(20) %>% 
  ggplot(aes(n, fct_reorder(go,n))) +
  geom_col(fill = 'darkblue') +
  # xlim(0, 50) +
  ggrepel::geom_label_repel(aes(label = description),
                            nudge_y = -0.5, nudge_x = 14, 
                            size = 3.5) +
  labs(y = 'GO terms',
       caption = 'Top 20 elements')

ggsave('../exploration/proteinfer/rare_prot_GO_filter.pdf',
       height = 6, width = 7)



# chromatine? -------------------------------------------------------------

full_blast %>% 
  filter()






