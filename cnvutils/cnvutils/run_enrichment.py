import altair as alt
import gseapy as gp
import gprofiler
import cptac
import cptac.pancan
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import warnings

def run_enrichr(input_file, cancer_type=None, gene_set_source="go"):
    """Run Enrichr. gene_set_source is either "go" or "reactome"."""

    input_df = pd.read_csv(input_file, sep="\t")
    
    if cancer_type is not None:
        cancer_df = input_df[input_df["cancers"].str.contains(cancer_type)]
        protein_list = cancer_df["protein"].tolist()
    else:
        protein_list = input_df["protein"].tolist()

    if gene_set_source == "go":
        gene_sets = ["GO_Biological_Process_2018"]

    elif gene_set_source == "reactome":
    
        path_here = os.path.abspath(os.path.dirname(__file__))
        gene_sets = os.path.join(path_here, "gene_set_libraries", "ReactomePathways.gmt")

    else:
        raise ValueError(f"Invalid value of '{gene_set_source}' for 'gene_set_source' parameter. Valid options are 'go' or 'reactome'.")

    enr = gp.enrichr(
        gene_list=protein_list,
        gene_sets=gene_sets,
        organism='Human',
        description='test_name',
        outdir=None,
        cutoff=0.05
    )
    
    return enr.res2d.sort_values(by="Adjusted P-value")

def run_gprofiler(input_file):

    input_df = pd.read_csv(input_file, sep="\t")
    protein_list = input_df["protein"].tolist()
    
    gp = gprofiler.GProfiler(return_dataframe=True)
    
    results = gp.profile(
        organism="hsapiens",
        query=protein_list,
        ordered=False,
        sources=["GO:BP", "KEGG", "REAC", "WP"]
    )
    
    return results
