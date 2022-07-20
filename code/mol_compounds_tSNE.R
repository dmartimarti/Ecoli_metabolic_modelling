library(tidyverse)

modelSeed_compounds = read_delim("tables/modelSeed_compounds.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)



polyamines_biocyc = read_delim("tables/polyamines_biocyc.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE) %>% 
  select(-Object_ID) 


modelSeed_compounds %>% 
  filter(inchikey %in% polyamines_biocyc$InChI_Key)



# step 1: add the polyamines compounds to the list

metformin = tibble(name = 'metformin', 
       smiles = 'CN(C)C(=N)N=C(N)N')

clean_modelSeed = modelSeed_compounds %>% 
  filter(!(inchikey %in% polyamines_biocyc$InChI_Key)) %>% 
  select(name, smiles) %>% 
  distinct(smiles, .keep_all = TRUE) %>% 
  bind_rows(polyamines_biocyc %>% 
              select(name, smiles)) %>% 
  bind_rows(metformin) %>% 
  drop_na()

clean_modelSeed %>% 
  rename(Drug = name) %>% 
  write.xlsx("tables/ModelSeed_polyamines_metf_smiles.xlsx")

# steps to do:
# add the remaining compounds to the complete list
# add metformin to the list
# use the pyMol2FuncGroup script to generate the list of func groups
# plot it with plotly, highlighting the polyamines and metformin 
