# All the functions

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

def load_params(path):
    """Load parameters for CNV analysis from the parameters file, created in the set_parameters_RUN_FIRST.ipynb notebook of the analysis folder for the current chromosome.

    Parameters:
    path (str): The path to the parameters file
    
    Returns:
    dict: The parameters
    """
    with open(path, "r") as file_obj:
        params = json.load(file_obj)

    return params

def load_tables(cancer_types, data_types, pancan):
    """Get the tables for the specified data types from the specified cancer types.

    Parameters:
    cancer_types (list of str): The cancer types to get data from
    data_types (list of str): The data types to get from each cancer, e.g. proteomics, CNV, transcriptomics, etc.
    pancan (bool): If False, use the regular cptac datasets. If true, use the cptac.pancan (harmonized) datasets.

    Returns:
    dict of str: dict of str: pd.DataFrame: A dict where the keys are data types and the values are dicts where the keys are cancer types and the values are dataframes of the proper data type.
    """
    # Initialize the dict we'll use to return tables
    all_tables = {}
    for data_type in data_types:
        all_tables[data_type] = {}

    # Load and save tables
    for cancer_type in cancer_types:
        cancer_type_tables = _load_cancer_type_tables(cancer_type, data_types, pancan)
        for data_type, df in cancer_type_tables.items():
            all_tables[data_type][cancer_type] = df

    return all_tables


def _load_cancer_type_tables(cancer_type, data_types, pancan):
    """Load the specified data tables from the given cancer type. We have this as a separate function instead of as part of load_tables so that the cancer dataset object will be allowed to be garbage collected after we're done with it, instead of sticking around and wasting RAM.

    Parameters:
    cancer_type (str): The cancer type to load
    data_types (list of str): The tables to get
    pancan (bool): If False, use the regular cptac datasets. If true, use the cptac.pancan (harmonized) datasets.

    Returns:
    dict of str: pd.DataFrame: The requested tables from the given cancer type, indexed by name.
    """

    # Load the cancer type
    if pancan:
        if cancer_type == "brca":
            ds = cptac.pancan.PancanBrca()
        elif cancer_type == "ccrcc":
            ds = cptac.pancan.PancanCcrcc()
        elif cancer_type == "colon":
            ds = cptac.pancan.PancanCoad()
        elif cancer_type == "endometrial":
            ds = cptac.pancan.PancanUcec()
        elif cancer_type == "gbm":
            ds = cptac.pancan.PancanGbm()
        elif cancer_type == "hnscc":
            ds = cptac.pancan.PancanHnscc()
        elif cancer_type == "lscc":
            ds = cptac.pancan.PancanLscc()
        elif cancer_type == "luad":
            ds = cptac.pancan.PancanLuad()
        elif cancer_type == "ovarian":
            ds = cptac.pancan.PancanOv()
        elif cancer_type == "pdac":
            ds = cptac.pancan.PancanPdac()
        else:
            raise ValueError(f"Invalid cancer type name '{cancer_type}'")

    else:
        if cancer_type == "brca":
            ds = cptac.Brca()
        elif cancer_type == "ccrcc":
            ds = cptac.Ccrcc()
        elif cancer_type == "colon":
            ds = cptac.Colon()
        elif cancer_type == "endometrial":
            ds = cptac.Endometrial()
        elif cancer_type == "gbm":
            ds = cptac.Gbm()
        elif cancer_type == "hnscc":
            ds = cptac.Hnscc()
        elif cancer_type == "lscc":
            ds = cptac.Lscc()
        elif cancer_type == "luad":
            ds = cptac.Luad()
        elif cancer_type == "ovarian":
            ds = cptac.Ovarian()
        elif cancer_type == "pdac":
            ds = cptac.Pdac()
        else:
            raise ValueError(f"Invalid cancer type name '{cancer_type}'")

    # Get the tables
    tables = {}

    for data_type in data_types:
        if pancan:
            if data_type == "CNV":
                tables[data_type] = ds.get_CNV()
            elif data_type == "proteomics":
                tables[data_type] = ds.get_proteomics(source="umich")
            elif data_type == "transcriptomics":
                tables[data_type] = ds.get_transcriptomics(source="bcm")
            else:
                raise ValueError(f"Invalid data type name '{data_type}'")
                
        else:
            tables[data_type] = ds._get_dataframe(data_type, tissue_type="both")

    return tables
    
def get_gene_locations():

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file = os.path.join(BASE_DIR, 'cnvutils', 'gene_locations.tsv')
    df = pd.read_csv(file, sep='\t', dtype={"chromosome": "O"}, index_col=[0,1], usecols=['Name', 'Database_ID', 'chromosome', 'start_bp', 'end_bp', 'arm'])

    return df

def get_event_genes(chrm, event_start, event_end, cis_or_trans):
    """Based on an event's start and end locations on a given chromosome, mark which genes are in the event and which ones aren't.

    Parameters:
    chrm (str): The chromosome the event is on.
    event_start (int): The base pair location of the event start.
    event_end (int): The base pair location of the event end.
    cis_or_trans (str): Either "cis" or "trans"; indicates whether you want proteins inside the event (cis) or outside the event (trans).

    Returns:
    pandas.Series: A boolean array indicating which genes are in the event and which ones aren't.
    """

    # Parameter processing
    cis_or_trans = cis_or_trans.lower()

    # Get locations
    locs = get_gene_locations()

    # Account for genes that go the opposite direction
    locs = locs.assign(
        first=locs[["start_bp", "end_bp"]].min(axis="columns"),
        last=locs[["start_bp", "end_bp"]].max(axis="columns")
    )

    # Create a filter for being in the event or not
    in_event = (
        (locs["chromosome"] == chrm) &
        (
            (
                (locs["first"] >= event_start) &
                (locs["first"] <= event_end)
            ) |
            (
                (locs["last"] >= event_start) &
                (locs["last"] <= event_end)
            )
        )
    )

    if cis_or_trans == "cis":
        ret = in_event[in_event]
    elif cis_or_trans == "trans":
        ret = in_event[~in_event]
    else:
        raise ValueError("Invalid value for 'cis_or_trans' parameter.")

    ret.name = "membership"

    ret = ret.\
    reset_index(drop=False).\
    drop(columns="membership").\
    drop_duplicates(keep="first").\
    sort_values(by=["Name", "Database_ID"])

    return ret

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

def make_cytoband_plot(chrm, show_xlabel=True, height=800):
    """Create a cytoband plot"""
    
    # Get cytoband info
    cytoband_data = get_cytoband_info()
    if chrm not in cytoband_data['#chromosome'].values:
        raise ValueError(f"Chromosome '{chrm}' not found in cytoband data. Make sure it's a string, not an int.")
    cytoband_data = cytoband_data[cytoband_data['#chromosome'] == chrm]
    
    # Create a column for colors
    cytoband_data = cytoband_data.assign(
        color=cytoband_data["stain"] + np.where(
            cytoband_data["density"].notna(),
            cytoband_data["density"].fillna(0).astype(int).astype(str),
            "",
        ),
    )

    # Make the chart
    bars = alt.Chart(cytoband_data).mark_bar().encode(
        x=alt.X(
            "#chromosome",
            title="Chromosome" if show_xlabel else None,
            axis=alt.Axis(
                labelAngle=0,
                labels=show_xlabel,
                ticks=show_xlabel,
                grid=False,
            ),
        ),
        y=alt.Y(
            "bp_start",
            stack=None,
            scale=alt.Scale(
                domain=[cytoband_data["bp_stop"].max(), 0],
                nice=False,
            ),
            title=None,
            axis=alt.Axis(
                # Hide the grid
                grid=False,
                
                # Set the tick frequency
                values=list(range(0, cytoband_data["bp_stop"].max(), 5000000)),
                
                # Below is a Vega expression that uses nested ternary if statements to return 0 if the
                # current tick mark's label is 0, otherwise check if it's a multiple of 10 million, and 
                # if so to trim off the first 8 characters (six zeros and two commas) and append "M" as 
                # an abbreviation for million base pairs, otherwise don't print a label for that tick mark
                labelExpr="datum.value == 0 ? 0 : datum.value % 10000000 ? null : slice(datum.label, -length(datum.label), -8) + 'Mb'",
            )
        ),
        color=alt.Color(
            "color",
            scale=alt.Scale(
                domain=["acen", "gneg",  "gpos25",    "gpos50",   "gpos75",  "gpos100", "gvar"],
                range=[ "crimson",  "white", "lightgray", "darkgray", "dimgray", "black",   "plum"],
            ),
            legend=None,
        ),
        order=alt.Order("bp_start", sort="ascending"),
    ).properties(height=height)
    
    return bars













def get_driver_genes():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file = os.path.join(BASE_DIR, 'cnvutils', 'bailey_driver_genes.csv')
    df = pd.read_csv(file, skiprows=3)
    return df


def get_cytoband_info():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file = os.path.join(BASE_DIR, 'cnvutils', 'NCBI_ideogram.csv')
    df = pd.read_csv(file)
    return df


def make_chromosome_plot(chromo, arm=None, start_bp=None, end_bp=None, genes=None, show_labels=True, title=None, ax=None, above=True):
    """ Create a cytoband plot and mark genes

    Parameters:
    chromo (str): The chromosome to be plotted
    arm (str): The chromosome arm to be plotted
    start_bp (int): the base pair to start plotting at
    end_bp (int): the base pair to end the plot
    genes (list or dict): a list of genes to mark on the plot; if using a dict, the key should be the color with the value as a list of genes to be marked in the given color.
    show_labels (bool): whether to show the gene names
    title (str): the title to show on the plot
    ax (Axes): the axes the plot should be generated on
    above (bool): If true labels will be placed above the plot. If false labels will be placed below.

    Results:
    Plot

    """
    cytoband_info = get_cytoband_info()
    data = cytoband_info[cytoband_info['#chromosome'] == chromo]
    locations = get_gene_locations()

    if arm:
        data = data[data.arm == arm]

    if start_bp:
        data = data[data.bp_stop > start_bp]
    else:
        start_bp = np.min(data.bp_start)

    if end_bp:
        data = data[data.bp_start < end_bp]
    else:
        end_bp = np.max(data.bp_stop)

    if above:
        label_location = 75
    else:
        label_location = 20

    colors = list()
    sections = list()
    for index, row in data.iterrows():
        sections.append((row['bp_start'], row['bp_stop']-row['bp_start']))
        if row['stain'] == 'gneg':
            colors.append('white')
        elif row['stain'] == 'gpos':
            if row['density'] == 25.0:
                colors.append('lightgray')
            elif row['density'] == 50.0:
                colors.append('gray')
            elif row['density'] == 75.0:
                colors.append('darkgray')
            else:
                colors.append('black')
        elif row['stain'] == 'acen':
            colors.append('red')
        else:
            colors.append('lightgray')
    if ax is None:
        fig, ax = plt.subplots()
        fig.set_figheight(0.5)
        fig.set_figwidth(30)
    ax.broken_barh(sections, (50,15), facecolors=colors, edgecolor="black")
    plt.axis('off')
    if title:
        ax.set_title(title, y=3.0, size='xx-large')
    not_found = list()
    if isinstance(genes, list):
        for gene in genes:

            loc = list(locations.loc[gene, 'start_bp'])[0]
            chromosome = list(locations.loc[gene, 'chromosome'])[0]
            if loc > start_bp and loc < end_bp and chromosome == chromo:
                ax.axvline(loc, 0, 15, color='r')
                if show_labels:
                    ax.text(loc, label_location, gene, rotation=90)

            else:
                not_found.append(gene)
    elif isinstance(genes, dict):
        for color in genes.keys():
            for gene in genes[color]:
                loc = list(locations.loc[gene, 'start_bp'])[0]
                chromosome = list(locations.loc[gene, 'chromosome'])[0]
                if loc > start_bp and loc < end_bp and chromosome == chromo:
                    ax.axvline(loc, 0, 15, color=color)
                    if show_labels:
                        ax.text(loc, label_location, gene, rotation=90)
                else:
                    not_found.append(gene)
    if len(not_found) > 0:
        warnings.warn(f'The following genes were not found within the event: {not_found}')

    return plt

def make_pvalue_plot(df, label_column, value_column, group_column=None, sort_column=None, sort_ascending=True, sig=0.05, show_sig=True, labels_per_plot=30):
    """
    @param df:
        The dataframe with pvalue information
    @param label_column:
        The name of the column that contains the labels for the x-axis of the figure
    @param value_column:
        The name of the column that contains the pvalues to be plotted
    @param group_column (optional):
        The name of the column that contains a grouping category. If provided, the groups will be indicated by color and
        a legend will be added to the figure
    @param sort_column (optional):
        The name of the column to sort the values on before plotting.
    @param sort_ascending
        Sort ascending vs. descending. This variable will only be used if sort_column is not None. Otherwise values will be
        plotted by position in the dataframe provided.
    @param sig:
        The significance value (before log transformation) to be plotted for comparison. The line can be turned off using the
        show_sig parameter.
    @param show_sig:
        Determines whether the significance line is shown in the plot.
    @param labels_per_plot:
        The number of labels on each plot. If there are more labels than fit on a single plot, additional plots will
        be created.
    """
    all_cancer_types = ['brca', 'colon', 'ccrcc', 'endo', 'gbm', 'hnscc', 'luad', 'lscc', 'ovarian']
    colors = sns.color_palette('tab10', n_colors=9)
    color_dict = dict()
    for i in range(9):
        color_dict[all_cancer_types[i]] = colors[i]
    df_copy = df.copy()
    df_copy['log_val'] = df_copy[value_column].apply(lambda x: -np.log10(x))
    if sort_column:
        df_copy = df_copy.sort_values(sort_column, ascending=sort_ascending)
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        x = list()
        for i in range(0, len(lst), n):
            x.append(lst[i:i + n])
        return x
    split_results = chunks(df_copy[label_column].unique(), labels_per_plot)
    def set_group(row):
        for i in range(len(split_results)):
            if row[label_column] in split_results[i]:
                return i
    df_copy['group_num'] = df_copy.apply(set_group, axis=1)
    if group_column:
        df_copy = df_copy.drop_duplicates(subset=[label_column, value_column, group_column])
    else:
        df_copy = df_copy.drop_duplicates(subset=[label_column, value_column])
#     palette = [color_dict[x] for x in df_copy.cancer.unique()]
    g = sns.FacetGrid(df_copy, row="group_num", aspect=4, sharex=False, sharey=False, legend_out=True)
    if group_column:
        plot = g.map_dataframe(sns.swarmplot, x=label_column, y="log_val", hue=group_column, palette=color_dict)
        g.add_legend()
    else:
        g.map_dataframe(sns.swarmplot, x=label_column, y="log_val", palette=color_dict)
    for ax in g.axes.ravel():
        ax.hlines(-np.log10(sig),*ax.get_xlim())
        ax.set_title("")
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.subplots_adjust(hspace=0.8, wspace=0.4)
    return g

def get_normal_expr_table():
    """Load the table of normal protein expression levels for different tissues. This table was downloaded from:
    https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20188503&file=msb188503-sup-0007-TableEV5.zip

    It was produced as part of this paper:
    Wang D, Eraslan B, Wieland T, et al. A deep proteome and transcriptome abundance atlas of 29 healthy human
    tissues. Mol Syst Biol. 2019;15(2):e8503. Published 2019 Feb 18. doi:10.15252/msb.20188503
    """
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(BASE_DIR, "cnvutils", "Table_EV5.xlsx")

    df = pd.read_excel(
        file_path,
        sheet_name="A. Protein copies"
    ).\
    rename(columns={
        "Gene name": "Gene_name",
        "Gene ID": "Gene_ID",
        "Protein ID": "Protein_ID"
    }).\
    set_index(["Gene_name", "Gene_ID", "Protein_ID"]).\
    sort_index().\
    reset_index(drop=False)

    return df
