# All the functions

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def get_gene_locations():

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file = os.path.join(BASE_DIR, 'cnvutils', 'gene_locations.tsv')
    df = pd.read_csv(file, sep='\t', dtype={"chromosome": "O"}, index_col=[0,1], usecols=['Name', 'Database_ID', 'chromosome', 'start_bp', 'end_bp', 'arm'])

    return df

def classify_event_genes(chrm, event_start, event_end, cis_or_trans):
    """Based on an event's start and end locations on a given chromosome, get either the cis or trans genes for the event.
    
    Parameters:
    chrm (str): The chromosome the event is on.
    event_start (int): The base pair location of the event start.
    event_end (int): The base pair location of the event end.
    cis_or_trans (str): "cis" to get the cis genes (genes contained in the event), or "trans" to get the trans genes (genes outside of the event)

    Returns:
    pandas.Series: A boolean array indicating which genes are in the event and which ones aren't.
    """
    
    # Parameter checking
    cis_or_trans = cis_or_trans.lower()
    if cis_or_trans != "cis" and cis_or_trans != "trans":
        raise ValueError(f"Valid options for cis_or_trans parameter are 'cis' or 'trans'. You passed '{cis_or_trans}'.")
    
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
    
    return in_event


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
