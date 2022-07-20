#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 12:15:26 2021

@author: dani
"""


from rdkit import Chem
import pandas as pd
from rdkit.Chem import Fragments
import re
from optparse import OptionParser 
from tqdm import tqdm

# create a OptionParser 
# class object 
parser = OptionParser() 

parser.add_option("-i", "--inp", 
                  dest = "input", 
                  help = "input file in form of an Excel file with two columns: 'Drugs' and 'smiles'",  
                  metavar = "INPUT")

parser.add_option("-o", "--out", 
                  dest = "output", 
                  help = "name of the output file that will store the chemical groups in a csv",  
                  metavar = "OUTPUT") 

(options, args) = parser.parse_args()


"""
It accepts a document with headers like:
Drug (can be any string)   smiles (smiles per compound)

It outputs a file with a matrix that has info about 
the functional groups present in each molecule
"""

def get_func_gr(molecule):
    '''
    Parameters
    ----------
    molecule : RDKit molecule type
        Molecule object from the RDKit library in python.

    Returns
    -------
    Returns a vector of the functional groups presence in the molecule.

    '''
    
    vector = []
    
    for func in frag_functions:
        str_to_eval = 'Fragments.'+func + '(molecule)'
        # print(str_to_eval)
        # print(func.replace('fr_',''))
        number = eval(str_to_eval)
        
        vector.append(number)
        
    return vector


### prep data table
drug_db = pd.read_excel(options.input)
drug_db  = drug_db.drop_duplicates('Drug') # remove duplicates
drug_db = drug_db[drug_db['smiles'].notna()] # remove na in smiles

# create a column with RDKit mols
drug_db['mol'] = drug_db['smiles'].map(Chem.MolFromSmiles)

drug_db = drug_db.set_index('Drug', drop = False)

# get the functions within Fragments and store them in a list
frag = dir(Fragments)
r = re.compile("fr_.*")
frag_functions = list(filter(r.match, frag)) # Read Note

# get column names
col_names = [i.replace('fr_','') for i in frag_functions]


def main():
    print(f'\n\nReading file {options.input} and saving the output in {options.output}.csv\n\n')

    drug_list = drug_db.Drug.to_list()

    drug_groups_total = []
    for drug in tqdm(drug_list, total = len(drug_list)):
        drug_groups = get_func_gr(drug_db.loc[drug].mol)
        drug_groups_total.append(drug_groups)

    global_dict = pd.DataFrame()

    df_drug = pd.DataFrame(drug_groups_total, index = drug_db.Drug, columns=col_names)

    df_drug.to_csv(f'{options.output}.csv', index=True)

    print('\n\nAnalysis done!\n\n')




if __name__ == "__main__":
    main()
