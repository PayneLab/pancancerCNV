import altair as alt
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import warnings

def make_ttest_plot(
    chromosome,
    arm,
    gain_or_loss,
    cis_or_trans,
    proteomics_or_transcriptomics,
    #cancer_types,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
):

    ttest_results = pd.\
    read_csv(os.path.join(
        data_dir,
        "data",
        f"chr{chromosome:0>2}{arm}_{gain_or_loss}_{cis_or_trans}_ttest.tsv"
    ), sep="\t").\
    rename(columns={"Name": "protein"})

    prots = ttest_results[ttest_results["adj_p"] <= 0.05].reset_index(drop=True)
    prots_cts = prots.groupby("cancer_type").count()[["protein"]]

    fail_prots = ttest_results[ttest_results["adj_p"] > 0.05].reset_index(drop=True)
    fail_cts = fail_prots.groupby("cancer_type").count()[["protein"]]

    prots_cts.insert(0, "count_type", "Significant difference")
    fail_cts.insert(0, "count_type", "No significant difference")

    counts = prots_cts.append(fail_cts).sort_index().reset_index(drop=False)

    event_effects_barchart = alt.Chart(counts).mark_bar().encode(
        x=alt.X(
            "count_type",
            axis=alt.Axis(
                title=None,
                labels=False,
                ticks=False,
            ),
            sort=["Significant difference"]
        ),
        y=alt.Y(
            "protein",
            axis=alt.Axis(
                title="Number of proteins"
            )
        ),
        color=alt.Color(
            "count_type",
            title=None,
            sort=["Significant difference"],
            scale=alt.Scale(
                domain=["Significant difference", "No significant difference"],
                range=["#2d3da4", "#d1d1d1"]
            )
        )
    ).facet(
        column=alt.Column(
            "cancer_type",
            title="Cancer type",
        )
    ).properties(
        title=f"Chr {chromosome}{arm} {cis_or_trans} {'protein' if proteomics_or_transcriptomics == 'proteomics' else 'RNA'} effects"
    ).configure_title(
        anchor="middle"
    ).configure_header(
        titleOrient="bottom",
        titlePadding=2,
        labelOrient="bottom",
    )

    return event_effects_barchart





































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
                # if so to trim off the last 8 characters (six zeros and two commas) and append "Mb" as 
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
