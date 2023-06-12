library(tidyverse)
library(readr)
library(cowplot)
library(here)
library(broom)
library(ggrepel)
library(glue)
library(openxlsx)
library(readxl)


theme_set(theme_cowplot(15))



# merge data from CLEAN ---------------------------------------------------


read_clean_splits = function(csv_file) {
  read_csv(csv_file, 
           col_names = c("gene", "predicted_EC")) %>% 
    separate_rows(predicted_EC, sep = ',') %>% 
    separate(predicted_EC, into = c("predicted_EC", "distance"), sep = '/') %>% 
    mutate(distance = as.numeric(distance)) 

}


read_clean_splits("raw_data/CLEAN/splits/split_ac_maxsep.csv")

root = "raw_data/CLEAN/splits/"

clean_ec = tibble()
for (csv in list.files(root, pattern = 'csv')) {
  temp_csv = read_clean_splits(paste0(root,csv))
  clean_ec = clean_ec %>% bind_rows(temp_csv)
}

# there are some rows that are extremely weird, let's fix the dataset
clean_ec %>% drop_na()


clean_extra = clean_ec %>% 
  separate_rows(X3, sep = ',') %>% 
  separate(X3, into = c("predicted_EC_X3", "distance_X3"), sep = '/') %>% 
  separate_rows(X4, sep = ',') %>% 
  separate(X4, into = c("predicted_EC_X4", "distance_X4"), sep = '/')

# get genes with x3 column info
x3_genes = clean_extra %>% 
  select(gene, predicted_EC_X3, distance_X3) %>% 
  drop_na() %>% 
  rename(predicted_EC = predicted_EC_X3,
         distance = distance_X3) %>% 
  mutate(distance = as.numeric(distance))

# get genes with x4 column info
x4_genes = clean_extra %>% 
  select(gene, predicted_EC_X4, distance_X4) %>% 
  drop_na() %>% 
  rename(predicted_EC = predicted_EC_X4,
         distance = distance_X4) %>% 
  mutate(distance = as.numeric(distance)) 


# get the final dataset
clean_ec = clean_ec %>% 
  select(gene, predicted_EC, distance) %>% 
  bind_rows(x3_genes, x4_genes) %>% 
  arrange(gene) %>% 
  mutate(predicted_EC = str_to_lower(predicted_EC))


enzyme_info = read_csv("raw_data/CLEAN/enzyme_information.csv") %>% 
  rename(predicted_EC = `EC Number`) %>% 
  mutate(predicted_EC = paste0('ec:', predicted_EC))


# exploration -------------------------------------------------------------

clean_ec = clean_ec %>% 
  left_join(enzyme_info)

clean_ec = clean_ec %>% rename(confidence = distance)


### IMPORTANT: ABOUT CONFIDENCE IN CLEAN TOOL
## https://clean.frontend.mmli1.ncsa.illinois.edu/configuration
## from running a few instances in their website, it's become clear to me
## that an acceptable confidence level starts from ~0.25, which leaves us 
## with only ~3500 genes with medium to high confidence levels. 
## all the other information must be used carefully 

clean_ec %>% filter(is.na(name))


clean_ec %>% 
  count(predicted_EC) %>% 
  arrange((n))


clean_ec %>% 
  count(gene) %>% 
  arrange(desc(n))

clean_ec %>%
  distinct(gene) %>% 
  count() 


clean_ec %>% 
  ggplot(aes(confidence)) +
  geom_histogram()

   
clean_ec %>% 
  filter(distance > 0.25) %>% 
  ggplot(aes(distance)) +
  geom_histogram()


   