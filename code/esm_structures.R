library(tidyverse)
library(NGLVieweR)
library(glue)
library(bio3d)


root = '/Users/danmarti/Documents/MRC_postdoc/My_projects/Ecoli_metabolic_modelling/prot_structures_esm/pdb_files'
proteins = list.files(root, pattern = 'pdb')

list.files('PG_protein_vis/data', pattern = 'pdb') %>% 
  write('PG_protein_vis/data/prot_list.txt')

prot_id = 3000
print(glue("Showing protein {proteins[prot_id]}"))
NGLVieweR(glue('{root}/{proteins[prot_id]}')) %>%
  addRepresentation("cartoon",
    param = list(name = "cartoon", colorScheme = "bfactor") 
  )

mol = read.pdb(glue("{root}/{proteins[prot_id]}"))

mol$atom %>% as_tibble() %>% 
  distinct(resid, resno) 


extract_sequence <- function(pdb_file) {
  pdb <- read.pdb(pdb_file)
  pdb = pdb$atom %>% as_tibble() %>% 
    distinct(resid, resno) 
}






