import cptac.utils
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

def select_genes_for_event(
    gene_locations,
    chromosome,
    event_start,
    event_end,
    cis_or_trans,
):

    # Take care of genes that go in the opposite direction of chromosome numbering
    gene_locations = gene_locations.assign(
        first=gene_locations[["start_bp", "end_bp"]].min(axis="columns"),
        last=gene_locations[["start_bp", "end_bp"]].max(axis="columns")
    )

    # Select genes on the chromosome
    chr_genes = gene_locations[gene_locations["chromosome"] == chromosome]

    # Create a filter for genes fully within the event
    in_event = ((chr_genes["first"] >= event_start) & (chr_genes["last"] <= event_end))

    # Select genes in the event for cis, or outside for trans
    if cis_or_trans == "cis":
        event_genes = chr_genes[in_event]
    elif cis_or_trans == "trans":
        event_genes = chr_genes[~in_event]
    else:
        raise ValueError("Invalid value for 'cis_or_trans' parameter.")

    # Drop our extra columns
    event_genes = event_genes.drop(columns=["first", "last"])

    return event_genes

def make_has_event_table(
    chromosome,
    arm,
    event_start,
    event_end,
    gain_or_loss,
    cancer_types,
    base_dir=os.getcwd(),
):
    # Get tables
    tables = load_input_tables(
        base_dir,
        data_types=["CNV"],
        cancer_types=cancer_types,
    )
    cnv_tables = tables["CNV"]
    
    # Get gene locations
    gene_locations = load_gene_locations(base_dir)

    for cancer_type in cnv_tables.keys():

        # Join in locations
        event_df = cnv_tables[cancer_type]
        event_df = event_df.transpose()
        if not isinstance(event_df.index, pd.MultiIndex):
            event_df = event_df.join(gene_locations.droplevel(1).drop_duplicates(keep="first"))
        else:
            event_df = event_df.join(gene_locations)
        
        # Slice out just the genes within the event
        event_df = select_genes_for_event(
            gene_locations=event_df,
            chromosome=chromosome,
            event_start=event_start,
            event_end=event_end,
            cis_or_trans="cis", # To get genes within the event
        )
        
        # Calculate gene lengths, drop other columns
        gene_lengths = event_df["end_bp"] - event_df["start_bp"]
        event_df = event_df.drop(columns=['chromosome', 'start_bp', 'end_bp', 'arm'])

        # Binarize all values to whether greater than/less than cutoff
        if gain_or_loss == "gain":
            bin_df = event_df.ge(INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).astype(int)
        elif gain_or_loss == "loss":
            bin_df = event_df.le(-INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).astype(int)
        else:
            raise ValueError(f"Invalid input '{gain_or_loss}' for gain_or_loss parameter")
        
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
            f"chr{chromosome:0>2}{arm}_{gain_or_loss}_{cancer_type}_has_event.tsv"
        ), sep='\t')

def event_effects_ttest(
    chromosome,
    arm,
    event_start,
    event_end,
    gain_or_loss,
    cis_or_trans,
    proteomics_or_transcriptomics,
    cancer_types,
    base_dir=os.getcwd(),
):
    # Get data_tables
    data_tables = load_input_tables(
        base_dir,
        data_types=[proteomics_or_transcriptomics],
        cancer_types=cancer_types,
    )
    data_tables = data_tables[proteomics_or_transcriptomics]
    
    # Get gene locations
    gene_locations = load_gene_locations(base_dir)

    # Select proteins
    # Note that the cnvutils.get_event_genes function uses Ensembl gene IDs for the 
    # Database_ID column, while the proteomics dataframes that have a Database_ID column
    # use RefSeq protein IDs. So, when we're selecting the genes we want, we ignore the 
    # Database_ID column if it is present, and just use gene names.
    selected_genes = select_genes_for_event(
        gene_locations=gene_locations,
        chromosome=chromosome,
        event_start=event_start,
        event_end=event_end,
        cis_or_trans=cis_or_trans,
    ).reset_index(drop=False)["Name"]

    for cancer_type in data_tables.keys():
        df = data_tables[cancer_type].transpose()
        
        if df.index.nlevels == 1:
            df = df[df.index.isin(selected_genes)]
        else:
            df = df[df.index.isin(selected_genes, level="Name")]

        data_tables[cancer_type] = df

    # Join in has_event data
    has_event = dict()
    for cancer_type in data_tables.keys():
        df = data_tables[cancer_type]
        df = df.transpose()
        event = pd.read_csv(os.path.join(
            base_dir,
            "data", 
            f"chr{chromosome:0>2}{arm}_{gain_or_loss}_{cancer_type}_has_event.tsv"
        ), sep='\t', index_col=0)
        event.index.rename('Name')
        df = df.join(event)
        df = df.dropna(subset=["event"])
        df = df.drop(columns="proportion")
        has_event[cancer_type] = df["event"]
        data_tables[cancer_type] = df
    
    # Run t-tests
    results_df = None
    for cancer_type in data_tables.keys():
        data_table = data_tables[cancer_type]
        results = cptac.utils.wrap_ttest(
            df=data_table, 
            label_column="event",
            correction_method="fdr_bh",
            return_all=True
        )
        results.set_index('Comparison', inplace=True)
        if isinstance(results.index[0], tuple):
            results[['Name', f'{cancer_type}_Database_ID']] = pd.DataFrame(
                results.index.values.tolist(),
                index=results.index
            )
            results.set_index(['Name', f'{cancer_type}_Database_ID'], inplace=True)
        else:
            results.index.name='Name'
        results.rename(columns={'P_Value': f'{cancer_type}_pvalue'}, inplace=True)
        if results_df is None:
            results_df = results
        else:
            results_df = results_df.join(results)
        
    # Append difference data
    diff_df = None
    for cancer_type in data_tables.keys():
        df = data_tables[cancer_type]
        df = df.drop("event", axis=1)
        results = df.apply(lambda x: _get_abundance_diff(x, has_event[cancer_type]))
        df = pd.DataFrame(results)
        if isinstance(df.index[0], tuple):
            df[['Name', f'{cancer_type}_Database_ID']] = pd.DataFrame(df.index.values.tolist(), index=df.index)
            df.set_index(['Name', f'{cancer_type}_Database_ID'], inplace=True)
        else:
            df.index.name='Name'
        df.rename(columns={0: f'{cancer_type}_diff'}, inplace=True)
        if diff_df is None:
            diff_df = df
        else:
            diff_df = diff_df.join(df)

    # Join the tables, reformat, and save
    results_df = results_df.join(diff_df)
    results_df = results_df.\
    reset_index(drop=False).\
    rename(columns={"Name": "protein"}).\
    set_index("protein")

    long_results = pd.DataFrame()

    for cancer_type in data_tables.keys():
        cancer_df = results_df.\
        loc[:, results_df.columns.str.startswith(cancer_type)].\
        dropna(axis="index", how="all").\
        reset_index(drop=False)

        # If the cancer type has database IDs, make a separate column that has them.
        # If not, create a column of NaNs (so that the tables all match)
        if f"{cancer_type}_Database_ID" in cancer_df.columns:
            cancer_df = cancer_df.rename(columns={f"{cancer_type}_Database_ID": "Database_ID"})
        else:
            cancer_df = cancer_df.assign(Database_ID=np.nan)

        # Rename the pvalue and diff columns to not have the cancer type
        cancer_df = cancer_df.rename(columns={
            f"{cancer_type}_pvalue": "adj_p",
            f"{cancer_type}_diff": "change"
        }).\
        assign(cancer_type=cancer_type)

        # Reorder the columns
        cancer_df = cancer_df[["cancer_type", "protein", "Database_ID", "adj_p", "change"]]

        # Append to the overall dataframe
        long_results = long_results.append(cancer_df)

    # Drop duplicate rows and reset the index
    long_results = long_results[~long_results.duplicated(keep=False)].\
    reset_index(drop=True)

    save_path = os.path.join(
        base_dir,
        "data", 
        f"chr{chromosome:0>2}{arm}_{gain_or_loss}_{cis_or_trans}_ttest.tsv"
    )
    long_results.to_csv(save_path, sep='\t', index=False)

# Helper functions
def _get_gain_counts(row):
    return (row > INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).sum()

def _get_loss_counts(row):
    return (row < -INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).sum()

def _get_abundance_diff(col, event):
    has_event = col[event]
    invert_list = [not x for x in event]
    no_event = col[invert_list]
    event_avg = has_event.mean()
    no_event_avg = no_event.mean()
    return event_avg - no_event_avg