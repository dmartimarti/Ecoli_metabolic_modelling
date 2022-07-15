
# libraries ---------------------------------------------------------------

library(tidyverse)
library(here)
library(cowplot)
library(sybil)
library(BacArena)

# sybil::SYBIL_SETTINGS("SOLVER","cplexAPI") # (optional)


# functions -------------------------------------------------------------


getMetaboliteProduction <- function(mod) {
  require(sybil)
  require(data.table)
  
  # MTF
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num])
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_cpd[0-9]+_e0", ex)])
  
  # FVA
  sol.fv <- fluxVar(mod, react = mod@react_id[grep("^EX_cpd[0-9]+_e0", mod@react_id)])
  
  dt <- data.table(ex       = rep(mod@react_id[grep("^EX_cpd[0-9]+_e0", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_cpd[0-9]+_e0", mod@react_id)],2),
                   dir      = c(rep("l",length(grep("^EX_cpd[0-9]+_e0", mod@react_id))),
                                rep("u",length(grep("^EX_cpd[0-9]+_e0", mod@react_id)))),
                   fv       = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(u>1e-6 & l >= 0)]
  
  dt <- merge(dt, dt.mtf, by = "ex")
  
  return(dt[order(-mtf.flux)])
}


# read metadata -------------------------


library(readxl)
metadata = read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx")

modelSeed_compounds <- read_delim("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/gapseq/modelSeed_compounds.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)


# get the files names
files_list = list.files(path = ".", pattern = ".RDS")

gnm_names = list.files(path = ".", pattern = ".RDS") %>% 
  str_sub(start = 1, end = -5)

# filter metadata for the bugs we want to analyse: AUS and ECOREF
meta_filt = metadata %>% 
  mutate(Genome = str_sub(fasta, start = 1, end = -7), .before = ID) %>% 
  filter(Discard == 'No') %>% 
  filter(Origin %in% c('AUS', 'ECOREF')) %>% 
  filter(Genome %in% gnm_names) %>% 
  distinct(Genome, .keep_all = T) %>% 
  mutate(theo_growth = 0) # make a new variable for the theoretical growth



# read models -------------------------------------------------------------


b1 <- readRDS("NT12060_237.RDS") # 
mg <- readRDS("NT12001_189.RDS") # for MG1655

getMetaboliteProduction(b1)[1:10]
getMetaboliteProduction(mg)[1:10]


b1_met = getMetaboliteProduction(b1)



# sybil sims --------------------------------------------------------------

### 1. read model ####

b1 <- readRDS("NT12060_237.RDS") 

### 2. run the FBA ####

optL = optimizeProb(b1, 
                    algorithm = 'fba', 
                    retOptSol= F)

### 3. get the obj function ####

optL$obj



# looping over the genomes ------------------------------------------------


growth = c()
for (model in files_list) {
  
  cat(glue::glue('Reading model {model}\n\n'))
  
  temp_model = readRDS(model)
  
  optL = optimizeProb(temp_model, 
                      algorithm = 'fba', 
                      retOptSol= F)
  
  growth = c(growth,optL$obj)
}

names(growth) = files_list

growth %>% 
  as_tibble(rownames = 'Model') %>% 
  ggplot(aes(x = fct_reorder(Model, value), y = value)) +
  geom_point()




# extract metabolites from models -----------------------------------------


tibble('met_name' = b1@met_name, 'met_id' = b1@met_id)

model_metabolites = tibble()
for (model in files_list) {
  
  cat(glue::glue('Reading model {model}\n\n'))
  
  temp_model = readRDS(model)
  
  temp_df = tibble('met_name' = temp_model@met_name, 
                   'met_id' = temp_model@met_id, 
                   'model' = model)
  
  model_metabolites = model_metabolites %>% bind_rows(temp_df)
}

# fix names and values


model_mets_fix = model_metabolites %>% 
  mutate(comp = str_sub(met_id, start = -3, end = -2),
         met_id = str_sub(met_id, start = 1, end = -5),
         model = str_sub(model, start = 1, end = -5)
    ) %>% 
  mutate(met_name = str_replace_all(met_name, '-c0', ''),
         met_name = str_replace_all(met_name, '-p0', ''),
         met_name = str_replace_all(met_name, '-e0', ''))

write_csv(model_mets_fix, '../model_metabolites.csv')

# model description -------------------------------------------------------

# media conditions
# Exchange reaction list

mg@met_id

# get cpd ids
id = str_sub(mg@met_id, end = -5)
compartments = str_sub(mg@met_id, start = -4)
# map to model seed database
new_cpd_id = tibble(id, compartments) %>% 
  inner_join(modelSeed_compounds %>% select(id, abbreviation)) %>% 
  mutate(new_vector = str_c(abbreviation, compartments)) %>% 
  pull(new_vector)
# change in the model
mg@met_id = new_cpd_id


# get exchange reactions

ex = findExchReact(mg)
opt = optimizeProb(mg)


# uptake react

upt = uptReact(ex)

ex[upt]


# genome coverage ---------------------------------------------------------

# let's try to calculate the genome coverage of these models

# start with an example
tibble(
  genes = b1@genes,
  rxn = b1@react_attr$rxn,
  start = b1@react_attr$sstart,
  end = b1@react_attr$send,
  status = b1@react_attr$status
) %>% 
  # drop_na(start, end) %>% 
  mutate(diff = end-start,
         direction = diff / abs(diff)) %>% view

## parallel version of data gathering ####
# parallel version of the loop, much faster
library(foreach)
library(doParallel)

# create the cluster
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# check if it is registered (optional)
foreach::getDoParRegistered()
print(my.cluster)


genes_size = foreach(i = 1:length(files_list),
        .packages=c("tidyverse"),
        .combine = 'rbind'
        ) %dopar% {

          # cat(glue::glue('Reading model {files_list[i]}\n\n'))

           temp_model = readRDS(files_list[i])

           temp_df = tibble(
             genes = temp_model@genes,
             rxn = temp_model@react_attr$rxn,
             start = temp_model@react_attr$sstart,
             end = temp_model@react_attr$send,
             status = temp_model@react_attr$status
           ) %>%
             drop_na(start, end) %>%
             mutate(diff = end-start,
                    direction = diff / abs(diff)) %>%
             mutate(Genome = files_list[i], .before = genes)

         }

stopCluster(my.cluster)


# 
# genes_size  = tibble()
# for (model in files_list) {
# 
#   cat(glue::glue('Reading model {model}\n\n'))
# 
#   temp_model = readRDS(model)
# 
#   temp_df = tibble(
#     genes = temp_model@genes,
#     rxn = temp_model@react_attr$rxn,
#     start = temp_model@react_attr$sstart,
#     end = temp_model@react_attr$send,
#     status = temp_model@react_attr$status
#   ) %>%
#     drop_na(start, end) %>%
#     mutate(diff = end-start,
#            direction = diff / abs(diff)) %>%
#     mutate(Genome = model, .before = genes)
# 
#   genes_size = genes_size %>% bind_rows(temp_df)
# }
# 
# 
# genes_size = genes_size %>% 
#   mutate(Genome = str_sub(Genome, start = 1, end = -5))


# A VERY IMPORTANT STEP
## 
genes_size_expand = genes_size %>% 
  # head(10000) %>%
  drop_na(diff) %>% 
  # filter(status != 'bad_blast') %>% 
  # select(Genome, genes, rxn, status) %>% 
  unnest(genes) %>% 
  distinct(Genome, genes, rxn, status, .keep_all = T) %>%
  # as genes are annotated at the end of the string, lets separate it
  separate(genes, into = c('new_start', 'new_end'), sep = ':',remove = F) %>% 
  # calculate new starts and ends. We'll need to fix the starts
  mutate(new_start = str_extract(new_start, '\\d+$'),
         new_start = as.numeric(new_start), # extracts the digits at the end of the string
         new_end = as.numeric(new_end)) %>% 
  mutate(new_diff = new_end - new_start,
         new_direction = new_diff / abs(new_diff)) %>% 
  mutate(new_diff = case_when(!is.na(new_diff) ~ new_diff,
                              is.na(new_diff) ~ diff),
         new_direction = case_when(!is.na(new_direction) ~ new_direction,
                                   is.na(new_direction) ~ direction)) %>% 
  select(Genome, genes, rxn,
         start = new_start, end = new_end, 
         diff = new_diff, direction = new_direction) 
  


# THERE ARE DUPLICATED REACTIONS! WTH!!
genes_size_expand %>% 
  # distinct(Genome, genes, rxn, .keep_all = T) %>%
  filter(rxn == 'ACETOOHBUTSYN-RXN') %>% 
  filter(Genome == '1.2')



genes_size_expand %>% 
  group_by(Genome) %>% 
  summarise(size = sum(abs(diff))) %>% 
  filter(size > 2.4e+06) %>% # there are a few outliers
  filter(size < 3.7e+06) %>% 
  ggplot(aes(x = size)) +
  geom_histogram(color = 'black',
                 fill = c("#2D51C4"),
                 alpha = 0.9) + 
  labs(
    x = 'Gene size sum',
    y = 'Count',
    caption = 'Values under 2e+06 and over 3.7e+06 have been filtered'
  ) +
  theme_cowplot(15)

ggsave('../exploration/genome_coverage.pdf',
       height = 4, width = 6)



## compare with real genome size ####
### absolute lengths ####

# Windows = "D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/assemblies/no_evo/quast_quality/transposed_report.tsv"
# Mac = "~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/assemblies/no_evo/quast_quality/transposed_report.tsv"
genome_info = read_delim("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/assemblies/no_evo/quast_quality/transposed_report.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         col_types = cols(Assembly = col_character()), 
                         trim_ws = TRUE) %>% 
  rename(genome = Assembly,
         genome_length = `Total length (>= 0 bp)`,
         genome_length_500 =`Total length`,
         gc_content = `GC (%)`) %>% 
  select(Genome = genome, genome_length, genome_length_500, gc_content)


genome_size_comp = genes_size_expand %>% 
  group_by(Genome) %>% 
  summarise(size = sum(abs(diff))) %>% 
  left_join(genome_info) %>% 
  # head(., 14) %>% 
  select(Genome,Metab_model_size=size, Genome_length= genome_length) %>% 
  mutate(Genome_diff = Genome_length - Metab_model_size) %>% 
  select(-Genome_length) %>% 
  pivot_longer(cols = Metab_model_size:Genome_diff,
               names_to = 'Fraction', 
               values_to = 'length')



genome_size_comp %>% 
  filter(!(Genome %in% c('6', 'SPC_3.1', '2.3','OP50'))) %>% 
  ggplot(aes(fill = Fraction, y = fct_reorder( Genome, length), x = length)) +
  geom_bar(position = 'stack', stat = 'identity', width = 1) +
  labs(x = 'Length',
       y = 'Genomes') +
  scale_fill_manual(name = NULL, 
                     # breaks = c(1,2),
                    values = c('#1BCBFA', '#FFA90F'),
                     labels = c("Total genome size",
                                "Genome covered")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank())


ggsave('../exploration/genome_coverage_genomeSize.pdf',
       height = 4, width = 6)


### percentage lengths ####

genome_size_comp_per = genes_size_expand %>% 
  group_by(Genome) %>% 
  summarise(size = sum(abs(diff))) %>% 
  left_join(genome_info) %>% 
  # head(., 14) %>% 
  select(Genome,Metab_model_size=size, Genome_length= genome_length) %>% 
  mutate(perc_covered = (Metab_model_size / Genome_length) * 100,
         perc_unc = 100 - perc_covered) %>% 
  select(-Genome_length, -Metab_model_size) %>% 
  arrange(desc(perc_covered)) %>% 
  mutate(Genome = factor(Genome, levels = Genome)) %>% 
  pivot_longer(cols = perc_covered:perc_unc,
               names_to = 'Fraction', 
               values_to = 'length') %>% 
  filter(!(Genome %in% c('6', 'SPC_3.1', '2.3','OP50'))) 


genome_size_comp_per %>% 
  mutate(new_label = case_when(Genome == 'NT12001_189' & 
                                 Fraction == 'perc_covered' ~ 'E. coli MG1655',
                               TRUE ~ '')) %>% 
  ggplot(aes(fill = factor(Fraction,
                           levels = c('perc_unc','perc_covered')), 
             y = Genome, x = length)) +
  geom_bar(position = 'stack', stat = 'identity', width = 1) +
  labs(x = '% covered by the metabolic model',
       y = 'Genomes') +
  ggrepel::geom_label_repel(aes(label = new_label),
                            max.overlaps = Inf,
                            box.padding = 0.4,
                            size = 3,
                            fill = 'white') +
  scale_fill_manual(name = NULL,
                    # breaks = c(1,2),
                    values = c('#1BCBFA', '#FFA90F'),
                    labels = c("Genome uncovered",
                               "Genome covered")) +
  scale_y_discrete(limits=rev) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank())

ggsave('../exploration/genome_coverage_genomeSize_percentage.pdf',
       height = 4, width = 6)








# BacArena ----------------------------------------------------------------



# Small fix to D/L-Lactate secretion (*) and model names
# bl <- rmReact(bl, react = "EX_cpd00221_e0")
b1@mod_desc <- "E. coli 12060"
mg@mod_desc <- "E. coli MG1655"

# Construct the organism objects for BacArena simulations
eure <- Bac(b1)
bilo <- Bac(mg)

# Construct the arena size 10x10 grid cells
arena <- Arena(n = 10, m = 10)

# For each organism, populate randomly 2 grid cells in the Arena as 
# 'starter culture'
arena <- addOrg(arena, eure, amount = 2)
arena <- addOrg(arena, bilo, amount = 2)

# add substrates to arena
arena_subs <- fread("../media/LBmed.csv") # same as gapfill medium
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]

arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
# Remove acetate from initial substrate list to see effect of Cross-Feeding
# arena <- rmSubs(arena, mediac = "EX_cpd00029_e0") 


# Simulation for 13 time steps
CF_sim <- simEnv(arena, time=10, sec_obj = "mtf", 
                 diff_par = T, cl_size = 6)

# Plot levels of Acetate, Buyrate, and Lactate as well as growth
par(mfrow=c(1,2))
plotCurves2(CF_sim,legendpos = "topleft",
            subs = c("cpd00211_e0","cpd00029_e0","cpd00159_e0", "cpd02799"),
            dict = list(cpd00011_e0 = "CO2", 
                        cpd00029_e0 = "Acetate", 
                        cpd00013_e0 = "Lactate",
                        cpd02799_e0 = 'Cosa'))


par(mfrow=c(1,2))
plotCurves2(CF_sim,legendpos = "topleft",
            subs = c("cpd02799_c0"),
            dict = list(cpd02799_c0 = 'Cosa'))


