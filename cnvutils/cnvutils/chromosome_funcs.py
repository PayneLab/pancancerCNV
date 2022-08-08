import cptac
import cptac.utils
import multiprocessing
import numpy as np
import os
import pandas as pd
import scipy
import warnings

from .constants import (
    ALL_CANCERS,
    INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF,
    PROPORTION_WITH_EVENT_CUTOFF,
    SIG_CUTOFF,
)
from .filenames import (
    get_cnv_counts_path,
    get_event_name,
    get_has_event_path,
    get_proportions_perm_test_results_path,
    get_ttest_results_path,
)
from .function_runners import multi_runner
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
    comparison,
    tissue_type=None,
    has_event=None,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
    save=True,
    permuted_event_data=None,
):

    if comparison not in ["tumor", "has_event"]:
        raise ValueError(f"Invalid value '{comparison}' for comparison parameter")

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

    to_drop = [] # For any empty dataframes
    for cancer_type in omics_dict.keys():
        
        df = omics_dict[cancer_type]

        # If applicable, select just the tissue type we want
        if comparison == "has_event":

            if tissue_type is None:
                raise ValueError(f"To compare based on has_event, you must choose a single tissue type.")
            elif tissue_type == "tumor":
                df = df[~df.index.str.endswith(".N")]
            elif tissue_type == "normal":
                df = df[df.index.str.endswith(".N")]
                df.index = df.index.str[:-2] # So that they'll match up
            else:
                raise ValueError(f"Invalid tissue_type '{tissue_type}'")

        if df.shape[0] == 0:
            to_drop.append(cancer_type)

        df = df.transpose().\
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

    # Drop any tables that didn't have any samples for the specified tissue type
    for cancer_type in to_drop:
        warnings.warn(f"Empty dataframe for cancer_type={cancer_type}, source={source}, level={level}, chromosome={chromosome}, arm={arm}, gain_or_loss={gain_or_loss}, cis_or_trans={cis_or_trans}, proteomics_or_transcriptomics={proteomics_or_transcriptomics}, comparison={comparison}, tissue_type={tissue_type}, has_event={has_event}.")
        del omics_dict[cancer_type]

    # Join in has_event data
    to_drop = [] # For any empty dataframes
    for cancer_type in omics_dict.keys():
        df = omics_dict[cancer_type]
        df = df.transpose()

        if permuted_event_data is None:
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
        else:
            event_name = get_event_name(
                source=source,
                level=level,
                chromosome=chromosome,
                arm=arm,
                gain_or_loss=gain_or_loss,
            )
            event = permuted_event_data[event_name][cancer_type]

        event.columns.name = "Name"
        event.columns = cptac.dataframe_tools.add_index_levels(
            to=event.columns,
            source=df.columns,
            fill="",
        )

        # If needed, create tissue_type column and get rid of suffixes on normal sample patient IDs
        if comparison == "tumor":

            # Create tissue_type column
            df = df.assign(tumor=~df.index.str.endswith(".N"))

            # Remove normal sample patient ID suffixes
            df.index = np.where(
                df.index.str.endswith(".N"),
                df.index.str[:-2],
                df.index,
            )

        # Join to event table
        df = df.\
        join(event, how="inner").\
        drop(columns=("proportion", ""))

        # If applicable, select just the event status that we want
        if comparison == "tumor":

            if has_event is None:
                raise ValueError(f"To compare based on tissue type, you must specify an event status with the 'has_event' parameter.")
            elif has_event:
                df = df[df["event"]]
            else:
                df = df[~df["event"]]

            # If there are no normal samples, we can't run a test
            if df["tumor"].value_counts().shape[0] == 1:
                to_drop.append(cancer_type)

            # Drop has_event column
            df = df.drop(columns=("event", ""))

        if df.shape[0] == 0:
            to_drop.append(cancer_type)

        omics_dict[cancer_type] = df

    # Drop any tables that didn't have any samples for the specified event status
    for cancer_type in to_drop:
        warnings.warn(f"Only one tissue type for cancer_type={cancer_type}, source={source}, level={level}, chromosome={chromosome}, arm={arm}, gain_or_loss={gain_or_loss}, cis_or_trans={cis_or_trans}, proteomics_or_transcriptomics={proteomics_or_transcriptomics}, comparison={comparison}, tissue_type={tissue_type}, has_event={has_event}.")
        del omics_dict[cancer_type]


    # Run t-tests
    if comparison == "tumor":
        label = "tumor"
    elif comparison == "has_event":
        label = "event"

    all_results = pd.DataFrame()
    for cancer_type in omics_dict.keys():

        omics = omics_dict[cancer_type]

        try:
            results = cptac.utils.wrap_ttest(
                df=omics, 
                label_column=omics[[label]].columns[0],
                alpha=SIG_CUTOFF,
                correction_method="fdr_bh",
                return_all=True,
                quiet=True,
            )

        except cptac.exceptions.InvalidParameterError as e:
            if str(e) == "No groups had enough members to pass mincount; no tests run.":
                warnings.warn(f"Too small sample size for cancer_type={cancer_type}, source={source}, level={level}, chromosome={chromosome}, arm={arm}, gain_or_loss={gain_or_loss}, cis_or_trans={cis_or_trans}, proteomics_or_transcriptomics={proteomics_or_transcriptomics}, comparison={comparison}, tissue_type={tissue_type}, has_event={has_event}.")
                continue
            else:
                raise e

        results = results.\
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
        [_get_has_event_sample_size, f"{comparison}_sample_size"],
        [_get_not_has_event_sample_size, f"not_{comparison}_sample_size"],
    ]:

        info_df = pd.DataFrame()
        for cancer_type in omics_dict.keys():

            df = omics_dict[cancer_type]
            label_col = df[label]
            df = df.drop(columns=label)

            results = df.apply(lambda x: func(x, label_col))

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

    if comparison == "has_event":
        comparison_name = "has_vs_not_has_event"
        group = tissue_type
    elif comparison == "tumor":
        comparison_name = "tumor_vs_normal"
        group = "has_event" if has_event else "not_has_event"

    save_path = get_ttest_results_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
        arm=arm,
        gain_or_loss=gain_or_loss,
        cis_or_trans=cis_or_trans,
        proteomics_or_transcriptomics=proteomics_or_transcriptomics,
        group=group,
        comparison_name=comparison_name,
    )

    if save:
        all_results.to_csv(save_path, sep='\t', index=False)
    else:
        return (save_path, all_results)

def get_has_vs_not_has_tumor_normal_diff_props(
    chromosomes_events,
    sources,
    levels,
    ttest_res=None,
    multi=True,
):

    # Get the proportions of genes affected and not affected for each event and category
    all_props = multi_runner(
        func=_get_ttest_sig_counts,
        sources=sources,
        levels=levels,
        chromosomes_events=chromosomes_events,
        more_dicts=[
            {
                "name": "ttest_res",
                "vals": [ttest_res],
            },
            {
                "name": "has_event",
                "vals": [True, False],
            },
            {
                "name": "proteomics_or_transcriptomics",
                "vals": ["proteomics", "transcriptomics"],
            },
            {
                "name": "cis_or_trans",
                "vals": ["cis", "trans"],
            },
        ],
        multi=multi,
    )

    all_props = pd.\
    concat(all_props, axis=0).\
    reset_index(drop=True)

    # Group proportions with the same event and categories, with has event paired with not has event
    all_props = all_props.\
    drop(columns="n").\
    infer_objects().\
    pivot(
        index=[
            "chromosome",
            "arm",
            "gain_or_loss",
            "cis_or_trans",
            "proteomics_or_transcriptomics",
            "cancer_types",
            "source",
            "level",
            "prop_name",
        ],
        columns="has_event",
        values="prop",
    ).\
    rename(columns={False: "not_has_prop", True: "has_prop"}).\
    reset_index(drop=False)

    all_props.columns.name = None

    # Split out the 3 types of proportions:
    #   - Proportion of genes with a tumor/normal p-value <= 0.05
    #   - Proportion of genes with a tumor/normal p-value > 0.05
    #   - Proportion of genes with a NaN tumor/normal p-value
    sig_props = all_props[all_props["prop_name"] == "sig_prop"]
    not_sig_props = all_props[all_props["prop_name"] == "not_sig_prop"]
    na_props = all_props[all_props["prop_name"] == "na_prop"]

    # Get a paired t-test p-value for all has event vs. not has event
    # proportions, then split into subgroups and recursively get p-values
    # for has vs. not has within those subgroups
    all_pvals = pd.DataFrame()
    for name, props in [
        ["sig_props", sig_props],
        ["not_sig_props", not_sig_props],
        ["na_props", na_props],
    ]:
        _, p_all = scipy.stats.ttest_rel(
            a=props["has_prop"],
            b=props["not_has_prop"],
            nan_policy="omit",
        )
        
        pvals = {
            "name": name,
            "p_all": p_all,
        }
        
        _split_groups(
            props, 
            split_cols=[
                "source",
                "proteomics_or_transcriptomics",
                "cis_or_trans",
                "gain_or_loss",
            ], 
            results=pvals,
            prefix="p"
        )
        
        pvals = pd.Series(pvals)
        all_pvals = pd.concat([all_pvals, pvals], axis=1)

    # Return all the p-values
    all_pvals = all_pvals.\
    transpose().\
    set_index("name", drop=True).\
    transpose().\
    reset_index(drop=False).\
    rename(columns={"index": "name"})

    all_pvals.columns.name = None

    # Reshape table to long format
    all_pvals = all_pvals.melt(
        id_vars=["name"],
        value_vars=["sig_props", "not_sig_props", "na_props"],
        var_name="group",
        value_name="adj_p",
    )

    return all_pvals

def permute_props(
    sources,
    levels,
    chromosomes_events,
    rng,
):
    # Optional fallback if we're not worried about reproducing exact permutations
    if rng is None:
        rng = np.random.default_rng()

    # Get permuted has_event tables for the events we're interested in
    all_event_perms = multi_runner(
        func=_permute_has_event,
        sources=["cptac", "gistic"],
        levels=["gene"],
        chromosomes_events={
            8: {
                "p": ["loss"],
                "q": ["gain"],
            },
        },
        more_dicts=[
            {
                "name": "rng",
                "vals": [rng],
            }
        ],
        multi=False,
    )
    
    # Convert tuples to dictionary
    all_event_perms = dict(all_event_perms)

    # Re-run all the tumor vs. normal t-tests, within samples with and without the event
    # separately, with the permuted has_event labels
    perm_res = multi_runner(
        func=event_effects_ttest,
        sources=["cptac", "gistic"],
        levels=["gene"],
        chromosomes_events={
            8: {
                "p": ["loss"],
                "q": ["gain"],
            },
        },
        more_dicts=[
            {
                "name": "save",
                "vals": [False],
            },
            {
                "name": "permuted_event_data",
                "vals": [all_event_perms],
            },
            {
                "name": "comparison",
                "vals": ["tumor"]
            },
            {
                "name": "has_event",
                "vals": [True, False],
            },
            {
                "name": "proteomics_or_transcriptomics",
                "vals": ["proteomics", "transcriptomics"],
            },
            {
                "name": "cis_or_trans",
                "vals": ["cis", "trans"],
            },
        ],
        multi=False,
    )

    # Convert tuple of (filename, df) into a dictionary
    # The next function will access the appropriate df based on what
    # its filename would have been if it was saved to disk
    perm_res = dict(perm_res)

    # Pass those t-test results to the proportion p value function and get the p values
    all_pvals = get_has_vs_not_has_tumor_normal_diff_props(
        chromosomes_events={
                8: {
                    "p": ["loss"],
                    "q": ["gain"],
                },
            },
        sources=["cptac", "gistic"],
        levels=["gene"],
        ttest_res=perm_res,
        multi=False,
    )

    # Return those p values so they can be added to the overall distribution
    return all_pvals

def props_permutation_test(
    n,
    sources,
    levels,
    chromosomes_events,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
):
    sq = np.random.SeedSequence()
    print(f"Entropy: '{sq.entropy}'")

    child_seeds = sq.spawn(n)
    args = [(
        sources,
        levels,
        chromosomes_events,
        np.random.default_rng(s),
    ) for s in child_seeds]
    
    with multiprocessing.Pool() as pool:
        results = pool.starmap(permute_props, args)
    
    all_pvals = pd.concat(results)

    save_path = get_proportions_perm_test_results_path(
        data_dir=data_dir,
        chromosome=f"{'_'.join([str(k) for k in chromosomes_events.keys()])}",
    )
    all_pvals.to_csv(save_path, sep="\t", index=False)

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

def _get_ttest_sig_counts(
    chromosome,
    arm,
    gain_or_loss,
    cis_or_trans,
    proteomics_or_transcriptomics,
    cancer_types,
    source,
    ttest_res,
    has_event=None,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
):
    
    res_path = get_ttest_results_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
        arm=arm,
        gain_or_loss=gain_or_loss,
        cis_or_trans=cis_or_trans,
        proteomics_or_transcriptomics=proteomics_or_transcriptomics,
        group="has_event" if has_event else "not_has_event",
        comparison_name="tumor_vs_normal",
    )

    if ttest_res is None:
        all_res = pd.read_csv(res_path, sep="\t")
    else:
        all_res = ttest_res[res_path]

    all_props = pd.DataFrame()
    for cancer_type in cancer_types:
        
        cancer_type_res = all_res[all_res["cancer_type"] == cancer_type]
        
        if cancer_type_res["adj_p"].size > 0:
            sig = (cancer_type_res["adj_p"] <= 0.05).sum()
            not_sig = (cancer_type_res["adj_p"] > 0.05).sum()
            na = cancer_type_res["adj_p"].isna().sum()

            sig_prop = sig / cancer_type_res["adj_p"].size
            not_sig_prop = not_sig / cancer_type_res["adj_p"].size
            na_prop = na / cancer_type_res["adj_p"].size
        else:
            sig_prop = np.nan
            not_sig_prop = np.nan
            na_prop = np.nan
            
        for prop in [
            ["sig_prop", sig_prop],
            ["not_sig_prop", not_sig_prop],
            ["na_prop", na_prop],
        ]:

            cancer_type_props = pd.Series({
                "chromosome": chromosome,
                "arm": arm,
                "gain_or_loss": gain_or_loss,
                "cis_or_trans": cis_or_trans,
                "proteomics_or_transcriptomics": proteomics_or_transcriptomics,
                "cancer_types": cancer_type,
                "source": source,
                "level": level,
                "has_event": has_event,
                "n": cancer_type_res["adj_p"].size,
                "prop_name": prop[0],
                "prop": prop[1],
            })

            all_props = pd.concat([all_props, cancer_type_props], axis=1)
        
    all_props = all_props.\
    transpose().\
    reset_index(drop=True)

    return all_props

def _split_groups(df, split_cols, results, prefix):
    """Recursively split proportions into subgroups and run t-tests on subgroups."""
    
    if len(split_cols) == 0:
        return
    
    col = split_cols[0]
    for val in df[col].unique():
        
        group = df[df[col] == val]
        
        _, p = scipy.stats.ttest_rel(
            a=group["has_prop"],
            b=group["not_has_prop"],
            nan_policy="omit",
        )
        
        new_prefix = "_".join([prefix, val])
        results[new_prefix] = p
        
        _split_groups(
            df=group,
            split_cols=split_cols[1:],
            results=results,
            prefix=new_prefix,
        )

def _permute_has_event(
    chromosome,
    arm,
    gain_or_loss,
    cancer_types,
    source,
    rng,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
):
    event_perms = {}
    
    for cancer_type in cancer_types:
        cancer_type_event_path = get_has_event_path(
            data_dir=data_dir,
            source=source,
            cancer_type=cancer_type,
            level=level,
            chromosome=chromosome,
            arm=arm,
            gain_or_loss=gain_or_loss,
        )

        cancer_type_event = pd.read_csv(cancer_type_event_path, sep='\t', index_col=0)
        perm_cancer_type_event = cancer_type_event.assign(event=rng.permutation(cancer_type_event["event"]))
        event_perms[cancer_type] = perm_cancer_type_event
        
    event_name = get_event_name(
        source=source,
        level=level,
        chromosome=chromosome,
        arm=arm,
        gain_or_loss=gain_or_loss,
    )

    return (event_name, event_perms)
