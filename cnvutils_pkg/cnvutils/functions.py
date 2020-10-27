# All the functions

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import warnings
import os

def get_gene_locations():

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file = os.path.join(BASE_DIR, 'cnvutils', 'gene_locations.tsv')
    df = pd.read_csv(file, sep='\t', dtype={"chromosome": "O"}, index_col=[0,1], usecols=['Name', 'Database_ID', 'chromosome', 'start_bp', 'end_bp', 'arm'])

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
    g = sns.FacetGrid(df_copy, row="group_num", aspect=4, sharex=False, sharey=False, legend_out=True)
    if group_column:
        g.map_dataframe(sns.swarmplot, x=label_column, y="log_val", hue=group_column, palette='muted')
        g.add_legend()
    else:
        g.map_dataframe(sns.swarmplot, x=label_column, y="log_val", palette="muted")
    for ax in g.axes.ravel():
        ax.hlines(-np.log10(sig),*ax.get_xlim())
        ax.set_title("")
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.subplots_adjust(hspace=0.8, wspace=0.4)
    return g



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

# def get_counts_table(cancer_type, dropna=True):
#     """
#     Returns the requested counts table from the data directory
#     """
#     counts_table = pd.read_csv(f'data/{cancer_type}_counts.csv', index_col=0)
#     if 'index' in counts_table.columns:
#         counts_table.set_index('index', inplace=True)
#         counts_table.index.name = 'Name'
#     else:
#         counts_table.set_index(['level_0', 'level_1'], inplace=True)
#         counts_table.index.rename(['Name', 'Database_ID'], inplace=True)
#     if dropna:
#         counts_table.dropna(subset=['chromo', 'location'], inplace=True)
#     return counts_table
#
# def get_ideogram_data(chromo):
#     """
#     Generates the information for an ideogram for the given chromosome
#     """
#     ideogram_data = pd.read_csv('data/NCBI_ideogram.csv') #cytoband info from NCBI
#     data = ideogram_data[ideogram_data['#chromosome'] == str(chromo)]
#     colors = list()
#     sections = list()
#     for index, row in data.iterrows():
#       sections.append((row['bp_start'], row['bp_stop']-row['bp_start']))
#       if row['stain'] == 'gneg':
#         colors.append('white')
#       elif row['stain'] == 'gpos':
#         if row['density'] == 25.0:
#           colors.append('lightgray')
#         elif row['density'] == 50.0:
#           colors.append('gray')
#         elif row['density'] == 75.0:
#           colors.append('darkgray')
#         else:
#           colors.append('black')
#       elif row['stain'] == 'acen':
#         colors.append('red')
#       else:
#         colors.append('lightgray')
#     return sections, colors
#
#
# def make_cnv_visual(cancer_types, chromo):
#     """
#     Makes a visual of the number of patients with gain and loss for the
#     types of cancer in cancer_types and the chromsome specified by chromo.
#     """
#     count_tables = list()
#     for cancer in cancer_types:
#         counts_table = get_counts_table(cancer)
#         counts_table = counts_table[counts_table.chromo == chromo]
#         counts_table['cancer'] = cancer
#         count_tables.append(pd.melt(counts_table, id_vars=['chromo', 'location', 'cancer'], value_vars=['gain', 'loss']))
#     count = 0
#     sections, colors = get_ideogram_data(chromo)
#     end_bp = sections[len(sections) - 1][0] + sections[len(sections) - 1][1]
#     fig, axs = plt.subplots(nrows=len(count_tables) + 1, sharex=True, sharey=True, num=0, figsize=(10,11), )
#     title = f'Chromosome {chromo}'
#     fig.suptitle(title, y=0.9, x=0.45)
#     plt.xlim(0,end_bp + (end_bp/5))
#     plt.ylim(0, 100)
#     for frame in count_tables:
#         axs[count].get_xaxis().set_visible(False)
#         axs[count].set_yticks([25,50])
#         axs[count].set_frame_on(False)
#         axs[count].text(end_bp + 5000000, 25, frame.cancer[0])
#         sns.lineplot(x="location", y="value", hue="variable", palette=['darkred', 'darkblue'], data=frame, ax=axs[count], legend=False)
#         if 'crossover' in frame.columns:
#             axs[count].axvline(frame.crossover[0], color='black', ls='dashed')
#         axs[count].set_ylabel("")
#         count += 1
#     plt.broken_barh(sections, (50,15), facecolors=colors, edgecolor="black")
#     plt.axis('off')
#     fig.legend(labels=("Amplifications", "Deletions"), loc='center right')
#     fig.text(0.07, 0.5, "Number of Samples", rotation="vertical")
