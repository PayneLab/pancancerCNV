ALL_CANCERS = [
    "brca",
    "ccrcc",
    "coad",
    "gbm",
    "hnscc",
    "lscc",
    "luad",
    "ov",
    "pdac",
    "ucec",
],

ALL_CHROMOSOMES = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "X",
    "Y",
]

# The CNV value magnitude required for a gene to count as amplified or deleted
INDIVIDUAL_GENE_CNV_MAGNITUDE_CUTOFF = 0.2

# The proportion of base pairs in an event region that must have an amplification for the sample to count
# as having an amplification event, or that must have a deletion for the sample to count as having a
# deletion event
PROPORTION_WITH_EVENT_CUTOFF = 0.8

# The cutoff below is the proportion of all patients in a cancer type that have 
# have a CNV event at a particular gene, for us to say that the gene is significantly
# gained or lost in that cancer type. This is used when defining event boundaries
# to decide which regions are gained or lost in each cancer type.
GENE_CNV_PROPORTION_CUTOFF =  0.2

# Parameters for saving charts
CHART_DPI = 600
CHART_FORMAT = "png"
CHART_RENDER_METHOD = "selenium" # npm doesn't allow for chart scaling
CHART_SCALE = 8.0

# Significance cutoff for B-H FDR adjusted p values
SIG_CUTOFF = 0.05
