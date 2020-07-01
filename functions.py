# All the functions
import pandas as pd

def get_counts_table(cancer_type, show_location=True):
    counts_table = pd.read_csv(f'data/{cancer_type}_counts.csv', index_col=0)
    return counts_table

def append_gene_locations():
    pass

def make_cnv_visual():
    pass

