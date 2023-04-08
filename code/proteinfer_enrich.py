# this script reads an input file (csv) with the following columns: gene, predicted_label, confidence, description, annotated
# it also reads a list of genes (one per line) from a csv file with a header "gene"
# the first csv file contains a list of genes and their GO terms annotations in the column "predicted_label"
# First filter the file to keep annotated == "Yes"
# Second, take the list of genes from the second file and calculate an enrichment test for the list of GO terms given in the first file
# The enrichment test is done using the hypergeometric test from scipy.stats.hypergeom
# The output is a csv file with the following columns: predicted_label, pvalue, enrichment
# use pandas to read the csv files
# import the files in the terminal line using the following command:
# python proteinfer_enrich.py input.csv genes.csv

import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from tqdm import tqdm
import multiprocessing as mp


# read the input file
df = pd.read_csv("/Users/danmarti/Documents/MRC_postdoc/My_projects/Ecoli_metabolic_modelling/tables/proteinfer_full.csv", sep=",", header=0)

# filter the input file to keep only annotated == "Yes"
df = df[df.annotated == "Yes"]

# read the list of genes
genes = pd.read_csv("/Users/danmarti/Documents/MRC_postdoc/My_projects/Ecoli_metabolic_modelling/tables/gene_result_list.csv", sep=",", header=0)

# get the list of GO terms
GO_terms = df["predicted_label"].unique()

# get the list of genes
genes_list = genes["gene"].unique()

# get the total number of genes
total_genes = len(genes['gene'].unique())

# get the total number of GO terms
total_GO_terms = len(GO_terms)

# get the total number of annotated genes
total_annotated_genes = len(df["gene"].unique())

# create a new dataframe to store the results
results = pd.DataFrame(columns=["predicted_label", "pvalue", "enrichment"])

# define a function to calculate the hypergeometric test
def hypergeom_test(GO_term):
    # get the list of genes annotated with the GO term
    GO_term_genes = df[df["predicted_label"] == GO_term]["gene"].unique()
    # get the number of genes annotated with the GO term
    n = len(GO_term_genes)
    # get the number of genes annotated with the GO term and in the list of genes
    k = len(np.intersect1d(GO_term_genes, genes_list))
    # calculate the hypergeometric test
    pvalue = hypergeom.sf(k-1, total_genes, n, total_annotated_genes)
    # calculate the enrichment
    enrichment = k/n
    # return the results
    return pd.Series({"predicted_label": GO_term, "pvalue": pvalue, "enrichment": enrichment})

# loop over the GO terms, use tqdm to show the progress, parallelize the loop using multiprocessing
# parallel for an M1 Macbook Pro
if __name__ == '__main__':
    with mp.Pool(processes=8) as pool:
        for result in tqdm(pool.imap_unordered(hypergeom_test, GO_terms), total=total_GO_terms):
            results = results.append(result, ignore_index=True)

# sort the results by pvalue
results = results.sort_values(by="pvalue")

# write the results to a csv file
results.to_csv("/Users/danmarti/Documents/MRC_postdoc/My_projects/Ecoli_metabolic_modelling/enrich_results.csv", sep=",", index=False)