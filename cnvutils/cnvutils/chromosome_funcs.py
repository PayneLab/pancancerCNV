import os
import pandas as pd

from .constants import ALL_CANCERS, INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF, PROPORTION_WITH_EVENT_CUTOFF
from .load_data import load_gene_locations, load_input_tables

def make_counts_table(
    chromosome,
    base_dir=os.getcwd(),
    cancer_types=ALL_CANCERS[0]
):
    # Get tables
    tables = load_input_tables(base_dir, data_types=["CNV"])
    cnv = tables["CNV"]
    
    # Get gene locations
    gene_locations = load_gene_locations(base_dir)
    chr_gene_locations = gene_locations[gene_locations["chromosome"] == chromosome]

    # Compile counts
    cnv_long = pd.DataFrame()
    for cancer_type in cancer_types:
        
        df = cnv[cancer_type].transpose()
        num_patients = df.shape[1]
        
        # Get just our chromosome
        df = df[df.index.get_level_values(0).isin(chr_gene_locations.index.get_level_values(0))]
        
        # Calculate counts
        df['gain'] = df.apply(_get_gain_counts, axis=1)
        df['loss'] = df.apply(_get_loss_counts, axis=1)
        
        # Join in locations
        df = df.join(chr_gene_locations)
        
        df = df.melt(
            id_vars=['start_bp', 'end_bp'], 
            value_vars=['gain', 'loss'], 
            ignore_index=False
        )
        
        df = df.assign(
            cancer_type_total_patients=num_patients,
            cancer=cancer_type
        )
        
        cnv_long = cnv_long.append(df)

    # Sort
    cnv_long = cnv_long.sort_values(['cancer', 'start_bp'])
    cnv_long = cnv_long.reset_index()

    # Save table
    cnv_long.to_csv(os.path.join(base_dir, "data", f"chr{chromosome:0>2}_cnv_counts.tsv"), sep='\t', index=False)

def make_has_event_table(
    chromosome,
    arm,
    event_start,
    event_end,
    loss_or_gain,
    cancer_types,
    base_dir=os.getcwd(),
):
    # Get tables
    tables = load_input_tables(base_dir, data_types=["CNV"])
    cnv_tables = tables["CNV"]
    
    # Get gene locations
    gene_locations = load_gene_locations(base_dir)

    for cancer_type in cnv_tables.keys():

        # Join in locations
        df = cnv_tables[cancer_type]
        df = df.transpose()
        if not isinstance(df.index, pd.MultiIndex):
            df = df.join(gene_locations.droplevel(1).drop_duplicates(keep="first"))
        else:
            df = df.join(gene_locations)

        # Get just our chromosome
        chr_df = df[df["chromosome"] == chromosome]
        
        # Slice out just the genes within the event
        event_df = chr_df[(chr_df.start_bp > event_start) & (chr_df.end_bp < event_end)]
        
        # Calculate gene lengths, drop other columns
        gene_lengths = event_df["end_bp"] - event_df["start_bp"]
        event_df = event_df.drop(columns=['chromosome', 'start_bp', 'end_bp', 'arm'])
        
        # Binarize all values to whether greater than/less than cutoff
        if loss_or_gain == "gain":
            bin_df = event_df.ge(INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).astype(int)
        elif loss_or_gain == "loss":
            bin_df = event_df.le(-INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).astype(int)
        else:
            raise ValueError(f"Invalid input '{loss_or_gain}' for loss_or_gain parameter")
        
        # Multiply every column by gene lengths
        scaled_df = bin_df.multiply(gene_lengths, axis="index")
        
        # Sum each column, and see what proportion of the overall event length it is. Keep if above cutoff.
        event_coding_length = gene_lengths.sum()
        proportions_df = scaled_df.sum(axis="index") / event_coding_length
        
        has_event = pd.DataFrame({
            "event": proportions_df >= PROPORTION_WITH_EVENT_CUTOFF,
            "proportion": proportions_df,
        })

        # Write to tsv
        has_event.to_csv(os.path.join(
            base_dir,
            "data", 
            f"chr{chromosome:0>2}{arm}_{loss_or_gain}_{cancer_type}_has_event.tsv"
        ), sep='\t')

# Helper functions
def _get_gain_counts(row):
    return (row > INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).sum()

def _get_loss_counts(row):
    return (row < -INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).sum()
