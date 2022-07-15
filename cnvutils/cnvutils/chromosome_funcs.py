import cptac
import cptac.utils
import numpy as np
import os
import pandas as pd

from .constants import (
    ALL_CANCERS,
    INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF,
    PROPORTION_WITH_EVENT_CUTOFF,
    SIG_CUTOFF,
)
from .filenames import (
    get_cnv_counts_path,
    get_has_event_path,
    get_ttest_results_path,
)
from .load_data import (
    get_cnv_counts,
    get_tables,
    get_ensembl_gene_locations,
    get_ncbi_gene_locations,
)

def make_counts_table(
    chromosome,
    source,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
    cancer_types=ALL_CANCERS[0],
):
    # Get tables
    if source == "cptac":

        id_name = "Database_ID"

        tables = get_tables(source=source, data_types=["CNV"], data_dir=data_dir)
        cnv = tables["CNV"]

        # Get gene locations
        gene_locations = get_ensembl_gene_locations(data_dir=data_dir)
        chr_gene_locations = gene_locations[gene_locations["chromosome"] == chromosome]

    elif source == "gistic":

        id_name = "NCBI_ID"

        tables = get_tables(source=source, data_types=[level], data_dir=data_dir)
        cnv = tables[level]

        if level == "gene":
            gene_locations = get_ncbi_gene_locations(data_dir=data_dir)
            chr_gene_locations = gene_locations[gene_locations["chromosome"] == chromosome]
        else:
            raise ValueError(f"Invalid level '{level}'")
    else:
        raise ValueError(f"Invalid data source '{source}'")
    
    # Compile counts
    cnv_long = pd.DataFrame()
    for cancer_type in cancer_types:
        
        df = cnv[cancer_type].transpose()
        num_patients = df.shape[1]

        # Get just our chromosome
        if source == "cptac" or source == "gistic" and level == "gene":
            df = df[df.index.get_level_values(id_name).isin(chr_gene_locations.index.get_level_values(id_name))]
        else:
            raise ValueError(f"Invalid level '{level}'")
        
        # Calculate counts
        df['gain'] = df.apply(_get_gain_counts, axis=1)
        df['loss'] = df.apply(_get_loss_counts, axis=1)
        
        # Join in locations
        if source == "cptac" or source == "gistic" and level == "gene":
            df = df.merge(
                right=chr_gene_locations,
                how="left",
                on=id_name,
            )
        else:
            raise ValueError(f"Invalid level '{level}'")
        
        df = df.melt(
            id_vars=['start_bp', 'end_bp'], 
            value_vars=['gain', 'loss'], 
            var_name="cnv_type",
            value_name="num_patients_with_cnv",
            ignore_index=False,
        )
        
        df = df.assign(
            cancer_type_total_patients=num_patients,
            cancer=cancer_type
        )

        cnv_long = pd.concat([cnv_long, df])

    # Sort
    cnv_long = cnv_long.sort_values(['cancer', 'start_bp'])
    cnv_long = cnv_long.reset_index()

    # Save table
    counts_path = get_cnv_counts_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
    )
    cnv_long.to_csv(counts_path, sep='\t', index=False)

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

    # Create a filter for genes fully within the event
    in_event = ((gene_locations["chromosome"] == chromosome) & (gene_locations["first"] >= event_start) & (gene_locations["last"] <= event_end))

    # Select genes in the event for cis, or outside for trans
    if cis_or_trans == "cis":
        event_genes = gene_locations[in_event]
    elif cis_or_trans == "trans":
        event_genes = gene_locations[~in_event]
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
    source,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
):
    # Get tables
    if source == "cptac":

        id_name = "Database_ID"

        tables = get_tables(source=source, data_types=["CNV"], data_dir=data_dir)
        cnv = tables["CNV"]

        # Get gene locations
        gene_locations = get_ensembl_gene_locations(data_dir=data_dir)

    elif source == "gistic":

        id_name = "NCBI_ID"

        tables = get_tables(source=source, data_types=[level], data_dir=data_dir)
        cnv = tables[level]

        if level == "gene":
            gene_locations = get_ncbi_gene_locations(data_dir=data_dir)
        else:
            raise ValueError(f"Invalid level '{level}'")
    else:
        raise ValueError(f"Invalid data source '{source}'")
    
    for cancer_type in cnv.keys():

        # Join in locations
        event_df = cnv[cancer_type]
        event_df = event_df.transpose()

        if source == "cptac" or source == "gistic" and level == "gene":
            event_df = event_df.merge(
                right=gene_locations,
                how="left",
                on=id_name,
            )
        else:
            raise ValueError(f"Invalid level '{level}'")
        
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
        event_df = event_df.drop(columns=[col for col in ['chromosome', 'start_bp', 'end_bp', 'arm', "Database_ID"] if col in event_df.columns])

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
        has_event_path = get_has_event_path(
            data_dir=data_dir,
            source=source,
            cancer_type=cancer_type,
            level=level,
            chromosome=chromosome,
            arm=arm,
            gain_or_loss=gain_or_loss,
        )
        has_event.to_csv(has_event_path, sep='\t')

def event_effects_ttest(
    chromosome,
    arm,
    event_start,
    event_end,
    gain_or_loss,
    cis_or_trans,
    proteomics_or_transcriptomics,
    cancer_types,
    source,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
):
    # Get omics tables
    omics_dict = get_tables(
        data_dir=data_dir,
        source="cptac",
        data_types=[proteomics_or_transcriptomics],
        cancer_types=cancer_types,
    )
    omics_dict = omics_dict[proteomics_or_transcriptomics]

    # Get gene locations
    gene_locations = get_ensembl_gene_locations(data_dir=data_dir)

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

    for cancer_type in omics_dict.keys():
        
        df = omics_dict[cancer_type].\
        transpose().\
        reset_index()

        df.columns.name = None

        # Get just the genes we want
        df = df[df["Name"].isin(selected_genes)]

        # Standardize index
        if "Database_ID" not in df.columns:
            df.insert(loc=1, column="Database_ID", value=np.nan)

        # Set the index again
        df = df.set_index(["Name", "Database_ID"])

        omics_dict[cancer_type] = df

    # Join in has_event data
    has_event = {}
    for cancer_type in omics_dict.keys():
        df = omics_dict[cancer_type]
        df = df.transpose()

        event_path = get_has_event_path(
            data_dir=data_dir,
            source=source,
            cancer_type=cancer_type,
            level=level,
            chromosome=chromosome,
            arm=arm,
            gain_or_loss=gain_or_loss,
        )
        event = pd.read_csv(event_path, sep='\t', index_col=0)

        event.columns.name = "Name"
        event.columns = cptac.dataframe_tools.add_index_levels(
            to=event.columns,
            source=df.columns,
            fill="",
        )

        df = df.join(event, how="inner")
        df = df.drop(columns="proportion")
        omics_dict[cancer_type] = df

        event = df[["event"]].copy(deep=True)
        event.columns = event.columns.get_level_values("Name")
        event = event["event"] # Make it a Series instead of a DataFrame

        has_event[cancer_type] = event

    # Run t-tests
    all_results = pd.DataFrame()
    for cancer_type in omics_dict.keys():

        omics = omics_dict[cancer_type]

        results = cptac.utils.\
        wrap_ttest(
            df=omics, 
            label_column=omics[["event"]].columns[0],
            alpha=SIG_CUTOFF,
            correction_method="fdr_bh",
            return_all=True,
            quiet=True,
        ).\
        set_index("Comparison").\
        rename(columns={"P_Value": "adj_p"}).\
        assign(cancer_type=cancer_type)

        results.index = pd.MultiIndex.from_tuples(
            tuples=results.index,
            names=["Name", "Database_ID"],
        )

        results = results.reset_index(drop=False)
        all_results = pd.concat([all_results, results])

    # Append difference data
    info_dfs = []
    for func, col_name in [
        [_get_abundance_change, "change"],
        [_get_has_event_sample_size, "has_event_sample_size"],
        [_get_not_has_event_sample_size, "not_has_event_sample_size"],
    ]:

        info_df = pd.DataFrame()
        for cancer_type in omics_dict.keys():
            df = omics_dict[cancer_type]
            df = df.drop("event", axis=1)

            results = df.apply(lambda x: func(x, has_event[cancer_type]))

            results = pd.DataFrame(results).\
            rename(columns={0: col_name}).\
            assign(cancer_type=cancer_type).\
            reset_index(drop=False)

            info_df = pd.concat([info_df, results])

        info_dfs.append(info_df)

    # Join the tables, reformat, and save
    for info_df in info_dfs:
        all_results = all_results.merge(
            right=info_df,
            how="outer",
            on=["Name", "Database_ID", "cancer_type"],
        )

    all_results = all_results.\
    reset_index(drop=True).\
    rename(columns={"Name": "protein"})

    save_path = get_ttest_results_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
        arm=arm,
        gain_or_loss=gain_or_loss,
        cis_or_trans=cis_or_trans,
        proteomics_or_transcriptomics=proteomics_or_transcriptomics,
    )
    all_results.to_csv(save_path, sep='\t', index=False)

# Helper functions
def _get_gain_counts(row):
    return (row > INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).sum()

def _get_loss_counts(row):
    return (row < -INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF).sum()

def _get_abundance_change(col, event):

    has_event = col[event]
    invert_list = [not x for x in event]
    no_event = col[invert_list]
    event_avg = has_event.mean()
    no_event_avg = no_event.mean()
    return event_avg - no_event_avg

def _get_has_event_sample_size(col, event):

    has_event = col[event]
    has_ct = has_event.notna().sum()

    return has_ct

def _get_not_has_event_sample_size(col, event):

    invert_list = [not x for x in event]
    no_event = col[invert_list]
    not_has_ct = no_event.notna().sum()

    return not_has_ct
