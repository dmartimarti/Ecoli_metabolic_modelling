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






prepareDatabase()

accessionToTaxa


