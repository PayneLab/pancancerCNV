import altair as alt
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import warnings

from .constants import (
    ALL_CANCERS,
    CHART_DPI,
    CHART_FORMAT,
    CHART_SCALE,
    GENE_CNV_PROPORTION_CUTOFF,
    SIG_CUTOFF,
)
from .filenames import (
    get_chr_gradient_plot_path,
    get_chr_line_plot_path,
    get_drivers_manhattan_plot_path,
    get_ttest_counts_plot_path,
    get_ttest_results_path,
)
from .load_data import (
    get_cnv_counts,
    get_cytoband_info,
    get_driver_genes,
    get_genes_ttest_results,
)

def make_chr_line_plot(
    source,
    chromosome,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
    cancer_types=ALL_CANCERS[0],
):

    # Turn off interactive plotting
    plt.ioff()

    # Get CNV counts
    cnv_counts = get_cnv_counts(
        source=source,
        level=level,
        chromosome=chromosome,
        data_dir=data_dir
    )

    # I was not able to find a good library for creating a visual of a chromosome
    # with their banding patterns, so I wrote this function to do it for me.

    ideogram_data = get_cytoband_info()
    chromo = ideogram_data[ideogram_data['chromosome'] == chromosome]
    colors = []
    sections = list()
    for index, row in chromo.iterrows():
        sections.append((row['bp_start'], row['bp_stop'] - row['bp_start']))
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

    # Set up our axes
    count = 0
    end_bp = sections[len(sections) - 1][0] + sections[len(sections) - 1][1]
    fig, axs = plt.subplots(nrows=len(cancer_types) + 1, sharex=True, sharey=True, num=0, figsize=(10,11))
    title = f"Chromosome {chromosome} CNV - {source}{' ' + level if level else ''}"
    fig.suptitle(title, y=0.9, x=0.45)
    plt.xlim(0,175138636)
    plt.xlim(0,end_bp + (end_bp/5))
    plt.ylim(0, 100)

    # Fill in the individual plots
    for cancer in cancer_types:
        frame = cnv_counts[cnv_counts["cancer"] == cancer]
        axs[count].get_xaxis().set_visible(False)
        axs[count].set_yticks([25,50])
        axs[count].set_frame_on(False)
        axs[count].text(end_bp + 5000000, 25, cancer)
        # Testing out the location event stuff
        axs[count].axvline(52110839, 0, 75, color='b')
        axs[count].axvline(202660, 0, 75, color='b')
        axs[count].axvline(37421341, 0, 75, color='b')
        sns.lineplot(
            x="start_bp", 
            y="num_patients_with_cnv", 
            hue="cnv_type", 
            palette=['darkred', 'darkblue'], 
            data=frame, 
            ax=axs[count], 
            legend=False)
        axs[count].set_ylabel("")
        count += 1
        
    # Add the chromosome banding map at the bottom
    plt.broken_barh(sections, (50,15), facecolors=colors, edgecolor="black")
    plt.axis('off')

    # Set up the legend
    red_line = mlines.Line2D([], [], color='darkred', label='Amplifications')
    blue_line = mlines.Line2D([], [], color='darkblue', label='Deletions')
    fig.legend(handles=[red_line, blue_line], loc='center right')

    # Set the Y axis label
    fig.text(0.07, 0.5, "Number of Samples", rotation="vertical")

    # Save the chart
    chart_path = get_chr_line_plot_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
        chart_format=CHART_FORMAT,
    )

    fig.savefig(chart_path, dpi=CHART_DPI, transparent=False, facecolor="white")
    plt.clf()

    return chart_path

def make_chr_gradient_plot(
    source,
    chromosome,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
    cancer_types=ALL_CANCERS[0],
):

    if source == "cptac":
        id_name = "Database_ID"
    elif source == "gistic":
        id_name = "NCBI_ID"
    else:
        raise ValueError(f"Invalid data source '{source}'")

    # Altair options
    alt.data_transformers.disable_max_rows()

    # Get CNV counts
    cnv_counts = get_cnv_counts(
        source=source,
        level=level,
        chromosome=chromosome,
        data_dir=data_dir
    )

    # Categorize counts
    gain_loss_counts = cnv_counts.pivot_table(index=[id_name, 'cancer'], columns='cnv_type')
    gain_loss_counts.columns = gain_loss_counts.columns.to_flat_index()
    gain_loss_counts = gain_loss_counts.drop(columns=[
        ('start_bp', 'gain'),
        ('end_bp', 'gain'),
        ('cancer_type_total_patients', 'gain'),
    ])
    gain_loss_counts = gain_loss_counts.rename(
        columns={
            ('cancer_type_total_patients', 'loss'): 'cancer_type_total_patients', 
            ('end_bp', 'loss'): 'end_bp', 
            ('start_bp', 'loss'): 'start_bp', 
            ('num_patients_with_cnv', 'gain'): 'gain', 
            ('num_patients_with_cnv', 'loss'): 'loss'}, 
    )

    gain_loss_counts['net_patient_ct'] = gain_loss_counts.gain - gain_loss_counts.loss

    gain_loss_counts = gain_loss_counts.\
    reset_index().\
    rename(columns={id_name: "gene"}).\
    sort_values(["cancer", "start_bp"])

    # Calculate a net_patient_prop column
    gain_loss_counts = gain_loss_counts.assign(
        net_patient_prop=gain_loss_counts["net_patient_ct"] / gain_loss_counts["cancer_type_total_patients"],
    )

    # Make the gradients
    grads = alt.Chart(gain_loss_counts).mark_rect().encode(
        x=alt.X(
            "cancer",
            title="Cancer type",
            axis=alt.Axis(
                labelAngle=40,
                ticks=False,
                grid=False,
                domain=False,
            ),
        ),
        y=alt.Y(
            "start_bp",
            title=None,
            axis=alt.Axis(
                labels=False,
                ticks=False,
                grid=False,
                domain=False,
                values=list(range(0, gain_loss_counts["end_bp"].astype(int).max(), 5000000)),
            )
        ),
        y2="end_bp",
        color=alt.Color(
            "net_patient_prop",
            title=[
                "Net proportion",
                "of patients with",
                "gain (red) or",
                "loss (blue)",
            ],
            scale=alt.Scale(
                scheme="redblue",
                domain=[-gain_loss_counts["net_patient_prop"].abs().max(), gain_loss_counts["net_patient_prop"].abs().max()],
                reverse=True,
            ),
        ),
    ).properties(
        width=300,
    )
    
    # Get the cytoband plot
    cytobands = make_cytoband_plot(
        chromosome,
        show_xlabel=False,
        width=20,
        height=500,
    )
    
    # Concatenate the plots
    grads = alt.hconcat(
        cytobands,
        grads,
        bounds="flush",
        title=f"Gene gain and loss on chromosome {chromosome}"
    ).resolve_scale(
        color="independent",
        y="shared",
    ).configure_title(
        anchor="middle"
    ).configure_view(
        stroke=None,
    )
    
    # Save the chart
    path = get_chr_gradient_plot_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
        chart_format=CHART_FORMAT,
    )

    grads.save(path, scale_factor=CHART_SCALE)

    return path

def find_gain_and_loss_regions(
    source,
    chromosome,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
    cancer_types=ALL_CANCERS[0],
):

    # Get our counts of number of patients with a CNV event at each gene
    counts = get_cnv_counts(
        source=source,
        level=level,
        chromosome=chromosome,
        data_dir=data_dir
    )

    # For each cancer type, calculate the minimum number of patients that need to have a CNV amplification
    # or deletion of a gene for us to consider that gene significantly amplified or deleted in that cancer
    # type, based on the GENE_CNV_PROPORTION_CUTOFF parameter. As of 12 Feb 2022, this is set at 20% of the total number
    # of patients in the cancer type.
    cutoffs = dict()

    for cancer_type in cancer_types:
        cutoffs[cancer_type] = counts[counts["cancer"] == cancer_type]["cancer_type_total_patients"].iloc[0] * GENE_CNV_PROPORTION_CUTOFF

    # Find gain and loss regions
    gains = _find_gain_or_loss_regions(
        counts=counts,
        cancer_types=cancer_types,
        gain_or_loss="gain",
        cutoffs=cutoffs,
    )

    losses = _find_gain_or_loss_regions(
        counts=counts,
        cancer_types=cancer_types,
        gain_or_loss="loss",
        cutoffs=cutoffs,
    )
    
    # Join the gain and loss data
    gains = gains.assign(event="gain")
    losses = losses.assign(event="loss", counts=losses["counts"] * -1)
    events = gains.append(losses)
    
    # Make the gains and losses plot
    events_plot = alt.Chart(events).mark_rect().encode(
        x=alt.X(
            "counts",
            title="Number of cancers with gain (positive) or loss (negative)",
        ),
        y=alt.Y(
            "start",
            title=None,
            axis=alt.Axis(
                labels=False,
                ticks=False,
                values=list(range(0, events["end"].astype(int).max(), 5000000)),
            )
        ),
        y2="end",
        color=alt.Color(
            "event",
            scale=alt.Scale(
                domain=["gain", "loss"],
                range=["darkred", "darkblue"],
            ),
        ),
        tooltip=alt.Tooltip(
            ["start", "end", "counts"],
            format=","
        ),
    )
    
    # Get the cytoband plot
    cytobands = make_cytoband_plot(
        chromosome,
        width=20,
        height=800,
    )
    
    # Concatenate the plots
    events_plot = alt.hconcat(
        cytobands,
        events_plot,
        bounds="flush"
    ).resolve_scale(
        color="independent",
        y="shared",
    )

    return events_plot, gains, losses

def make_ttest_counts_plot(
    chromosome,
    arm,
    gain_or_loss,
    cis_or_trans,
    proteomics_or_transcriptomics,
    source,
    level=None,
    data_dir=os.path.join(os.getcwd(), "..", "data"),
):

    # Get our ttest results file
    ttest_results_path = get_ttest_results_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
        arm=arm,
        gain_or_loss=gain_or_loss,
        cis_or_trans=cis_or_trans,
        proteomics_or_transcriptomics=proteomics_or_transcriptomics,
    )

    ttest_results = pd.\
    read_csv(ttest_results_path, sep="\t").\
    rename(columns={"Name": "protein"})

    prots = ttest_results[ttest_results["adj_p"] <= 0.05].reset_index(drop=True)
    prots_cts = prots.groupby("cancer_type").count()[["protein"]]

    fail_prots = ttest_results[ttest_results["adj_p"] > 0.05].reset_index(drop=True)
    fail_cts = fail_prots.groupby("cancer_type").count()[["protein"]]

    prots_cts.insert(0, "count_type", "Significant difference")
    fail_cts.insert(0, "count_type", "No significant difference")

    counts = pd.concat([prots_cts, fail_cts]).sort_index().reset_index(drop=False)

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
        title=f"{source} {level + ' level ' if level else ''}chr {chromosome}{arm} {cis_or_trans} {'protein' if proteomics_or_transcriptomics == 'proteomics' else 'RNA'} effects"
    ).configure_title(
        anchor="middle"
    ).configure_header(
        titleOrient="bottom",
        titlePadding=2,
        labelOrient="bottom",
    )

    # Save the chart
    path = get_ttest_counts_plot_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
        arm=arm,
        gain_or_loss=gain_or_loss,
        cis_or_trans=cis_or_trans,
        proteomics_or_transcriptomics=proteomics_or_transcriptomics,
        chart_format=CHART_FORMAT,
    )

    event_effects_barchart.save(path, scale_factor=CHART_SCALE)

    return path

def make_genes_manhattan_plot(
        genes,
        source,
        chromosome,
        arm,
        gain_or_loss,
        cis_or_trans,
        proteomics_or_transcriptomics,
        level=None,
        facet_row=None,
        data_dir="../data",
        title=None,
):

    ttest_results = get_genes_ttest_results(
        genes=genes,
        source=source,
        chromosome=chromosome,
        arm=arm,
        gain_or_loss=gain_or_loss,
        cis_or_trans=cis_or_trans,
        proteomics_or_transcriptomics=proteomics_or_transcriptomics,
        level=level,
        data_dir=data_dir,
    )

    ttest_results = ttest_results.assign(
        line=-np.log10(SIG_CUTOFF),
        neg_log_adj_p=-np.log10(ttest_results["adj_p"]),
    )

    base = alt.Chart(ttest_results)

    dots = base.mark_point().encode(
        x=alt.X(
            "cancer_type",
            title=None,
        ),
        y=alt.Y(
            "neg_log_adj_p",
            title="-log(p)",
        ),
        color=alt.Color(
            "cancer_type",
            title="Cancer type",
        ),
    )

    line = base.mark_rule(color="crimson").encode(
        y="line",
    )

    if title is None:
        title = f"{source} {level + ' level ' if level else ''}chr {chromosome}{arm} {cis_or_trans} {'protein' if proteomics_or_transcriptomics == 'proteomics' else 'RNA'} effects"

    mplot = alt.layer(
        dots,
        line,
    ).facet(
        facet="protein",
        columns=7,
        spacing=alt.RowColnumber(column=0, row=20),
    ).properties(
        title=title,
    ).configure_title(
        anchor="middle"
    ).configure_header(
        title=None,
    )

    return mplot

def make_drivers_manhattan_plot(
    source,
    chromosome,
    arm,
    gain_or_loss,
    cis_or_trans,
    proteomics_or_transcriptomics,
    level=None,
    data_dir="../data",
):

    if source == "cptac":
        locs_source = "ensembl"
    elif source == "gistic":
        locs_source = "ncbi"
    else:
        raise ValueError(f"Invalid 'source' parameter '{source}'")
    
    genes = get_driver_genes(locs_source=locs_source)

    mplot = make_genes_manhattan_plot(
            genes=genes["Name"],
            source=source,
            chromosome=chromosome,
            arm=arm,
            gain_or_loss=gain_or_loss,
            cis_or_trans=cis_or_trans,
            proteomics_or_transcriptomics=proteomics_or_transcriptomics,
            level=level,
            data_dir=data_dir,
            title=f"{source} {level + ' level ' if level else ''}chr {chromosome}{arm} {cis_or_trans} {'protein' if proteomics_or_transcriptomics == 'proteomics' else 'RNA'} drivers effects"
    )

    # Save the chart
    path = get_drivers_manhattan_plot_path(
        data_dir=data_dir,
        source=source,
        level=level,
        chromosome=chromosome,
        arm=arm,
        gain_or_loss=gain_or_loss,
        cis_or_trans=cis_or_trans,
        proteomics_or_transcriptomics=proteomics_or_transcriptomics,
        chart_format=CHART_FORMAT,
    )

    mplot.save(path, scale_factor=CHART_SCALE)

    return path

def make_cytoband_plot(chrm, width, height, show_xlabel=True):
    """Create a cytoband plot"""
    
    # Get cytoband info
    cytoband_data = get_cytoband_info()
    cytoband_data = cytoband_data[cytoband_data['chromosome'] == chrm]
    
    # Create a column for colors
    cytoband_data = cytoband_data.assign(
        color=cytoband_data["stain"] + np.where(
            cytoband_data["density"].notna(),
            cytoband_data["density"].fillna(0).astype(int).astype(str),
            "",
        ),
    )

    # Make the chart
    bars = alt.Chart(cytoband_data).mark_rect().encode(
        x=alt.X(
            "chromosome",
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
        y2="bp_end",
        color=alt.Color(
            "color",
            scale=alt.Scale(
                domain=["acen", "gneg",  "gpos25",    "gpos50",   "gpos75",  "gpos100", "gvar"],
                range=[ "crimson",  "white", "lightgray", "darkgray", "dimgray", "black",   "plum"],
            ),
            legend=None,
        ),
        order=alt.Order("bp_start", sort="ascending"),
    ).properties(
        width=width,
        height=height,
    )
    
    return bars

def _find_gain_or_loss_regions(
    counts,
    cancer_types,
    gain_or_loss,
    cutoffs,
):
    event_locations = dict()
    for cancer in cancer_types:
        
        df = counts[(counts.cnv_type == gain_or_loss) & (counts.cancer == cancer)].sort_values('start_bp')
        values = list(df.num_patients_with_cnv)
        events = list()
        start = None
        for i in range(0, len(values)):
            val = values[i]
            if val > cutoffs[cancer]:
                if start is None:
                    start = i
            else:
                if start is not None:
                    events.append((start, i))
                    start = None
        if start is not None:
            events.append((start, len(values)-1))
        event_locations[cancer] = []
        for event in events:
            start_bp = df.iloc[event[0]].start_bp
            end_bp = df.iloc[event[1]].start_bp
            event_locations[cancer].append((start_bp, end_bp - start_bp))

    event_patients = list()
    for cancer in event_locations.keys():
        events = event_locations[cancer]
        for event in events:
            start = event[0]
            end = event[0] + event[1]
            event_patients.append((start, 1))
            event_patients.append((end, 0))
    event_patients.sort()

    count = 0
    current_bp = 0
    start = list()
    end = list()
    size = list()
    total = list()
    for patient in event_patients:
        if patient[0] != current_bp:
            start.append(current_bp)
            end.append(patient[0])
            size.append(patient[0] - current_bp)
            total.append(count)
            current_bp = patient[0]
        if patient[1] == 1:
            count += 1
        else:
            count -= 1
    event_data = pd.DataFrame({'start': start, 'end': end, 'counts': total, 'length': size}).sort_values('start')

    return event_data

