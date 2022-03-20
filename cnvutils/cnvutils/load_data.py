import cptac
import cptac.pancan
import glob
import json
import numpy as np
import os
import pandas as pd
import pyensembl
import requests
import sys

from .constants import ALL_CANCERS

def save_input_tables(pancan, base_dir=os.getcwd()):
    """Load CNV, transcriptomics, and proteomics tables for all cancers, create
    gene locations table, and save all of them.
    """

    # Load tables
    tables = _load_cptac_tables(
        cancer_types=ALL_CANCERS[0],
        data_types=[
            "CNV",
            "transcriptomics",
            "proteomics",
        ],
        pancan=pancan,
    )

    # Get list of all genes we have CNV data for. These are the ones we'll need the locations of.
    genes = None
    for df in tables["CNV"].values():
        if genes is None:
            genes = df.columns.copy(deep=True) # Just for the first one
        else:
            genes = genes.union(df.columns)

    genes = genes.to_frame()
    gene_locations, not_found_genes = query_gene_locations_database(genes)

    # Create a data directory in the directory the function was called from
    input_data_dir = os.path.join(base_dir, "data", "sources")
    cptac_data_dir = os.path.join(input_data_dir, "cptac_tables")
    os.makedirs(input_data_dir, exist_ok=True)
    os.makedirs(cptac_data_dir, exist_ok=True)

    # Save the gene locations
    gene_locations_save_path = os.path.join(input_data_dir, "gene_locations.tsv.gz")
    gene_locations.to_csv(gene_locations_save_path, sep="\t")

    # Save the omics tables
    for data_type, cancer_types in tables.items():
        for cancer_type, df in cancer_types.items():

            file_name = f"{cancer_type}_{data_type}.tsv.gz"
            save_path = os.path.join(cptac_data_dir, file_name)
            df.to_csv(save_path, sep="\t")

def load_input_tables(base_dir, data_types=["CNV", "proteomics", "transcriptomics"], cancer_types=ALL_CANCERS[0]):

    # Standardize data_types and cancer_types
    data_types = [data_type.lower() for data_type in data_types]
    cancer_types = [cancer_type.lower() for cancer_type in cancer_types]

    # Get the data tables directory
    cptac_tables_dir = os.path.join(base_dir, "data", "sources", "cptac_tables")

    # Get list of tables to load
    all_table_paths = sorted(glob.glob(os.path.join(cptac_tables_dir, "*")))
    load_table_paths = []
    for path in all_table_paths:
        cancer_type, data_type = path.split(os.sep)[-1].split(".")[0].split("_")
        if cancer_type.lower() in cancer_types and data_type.lower() in data_types:
            load_table_paths.append(path)

    # Load the tables
    tables = {}
    for i, path in enumerate(load_table_paths):
        cancer_type, data_type = path.split(os.sep)[-1].split(".")[0].split("_")

        print(f"Loading {cancer_type} {data_type} ({i + 1}/{len(load_table_paths)})...{' ' * 30}", end="\r")

        if data_type not in tables.keys():
            tables[data_type] = {}

        tables[data_type][cancer_type] = pd.read_csv(
            path,
            sep="\t",
            index_col=0,
            header=[0, 1]
        )

    # Clear last loading message
    print(" " * 80, end="\r")

    return tables

def load_gene_locations(base_dir=os.getcwd()):

    gene_locations_path = os.path.join(base_dir, "data", "sources", "gene_locations.tsv.gz")
    gene_locations = pd.read_csv(gene_locations_path, sep="\t", index_col=[0, 1])

    return gene_locations

def query_gene_locations_database(genes):

    # Find latest Ensembl release supported by pyensembl
    latest_release = 100
    invalid_release = False
    while not invalid_release:
        try:
            pyensembl.ensembl_release_versions.check_release_number(latest_release + 1)
        except ValueError:
            invalid_release = True
        else:
            latest_release += 1

    # Load pyensembl Ensembl API
    ensembl = load_ensembl_release(latest_release)

    gene_names = []
    db_id = []
    chromosome = []
    start_bp = []
    end_bp = []

    # Look up genes using the latest Ensembl release
    db_not_found = []
    no_db_id = []
    name_not_found = []
    for gene in genes.itertuples(index=False):

        info = None # Reset from previous iteration

        if gene.Database_ID:
            try:
                info = ensembl.gene_by_id(gene.Database_ID.split(".")[0])
            except ValueError:
                db_not_found.append(gene)
        else:
            no_db_id.append(gene)
            try:
                info = ensembl.genes_by_name(gene.Name)[0]
            except ValueError:
                name_not_found.append(gene)

        if info is not None: # If info is not None, then the lookup was successful
            gene_names.append(gene.Name)
            db_id.append(gene.Database_ID)
            chromosome.append(info.contig)
            start_bp.append(info.start)
            end_bp.append(info.end)

    # For genes with database IDs who weren't found in the latest Ensembl release, we can query the Ensembl API to see if they're from an older Ensembl release, and then look them up with those older releases
    # First we'll look up which release each identifier is from
    old_db_release_not_found = []
    old_db_genes = []
    old_db_release_nums = []
    for i, gene in enumerate(db_not_found):

        print(f"Looking up release for ID {i + 1}/{len(db_not_found)}", end="\r")
        try:
            old_release_number = _lookup_old_release(gene.Database_ID.split(".")[0])
        except ValueError:
            old_db_release_not_found.append(gene)
        else:
            old_db_genes.append(gene)
            old_db_release_nums.append(old_release_number)

    # Second, for the IDs we did find an old release for, we'll sort them by release and then look up all the IDs for each release
    old_db_df = pd.\
    DataFrame({
        "gene": old_db_genes,
        "old_release_number": old_db_release_nums,
    })

    old_db_dict = dict(tuple(old_db_df.groupby("old_release_number")["gene"]))

    old_db_not_found = []
    release_too_new = []
    for i, (old_release_number, old_genes) in enumerate(old_db_dict.items()):

        if old_release_number > latest_release:
            release_too_new.extend(old_genes.tolist())
            continue

        old_ensembl = load_ensembl_release(old_release_number)

        print(f"Looking up {len(old_genes)} genes with release {old_release_number} (release {i + 1}/{len(old_db_dict.keys())})     ", end="\r")
            
        for j, gene in enumerate(old_genes):

            try:
                info = old_ensembl.gene_by_id(gene.Database_ID.split(".")[0])
            except ValueError as e:
                old_db_not_found.append(gene)
            else:
                gene_names.append(gene.Name)
                db_id.append(gene.Database_ID)
                chromosome.append(info.contig)
                start_bp.append(info.start)
                end_bp.append(info.end)

    # Clear message
    print(" " * 80, end="\r")

    # If we couldn't find the database ID in an older version, we'll last try looking it up by name
    name_and_db_not_found = []
    for gene in old_db_release_not_found + old_db_not_found:

        try:
            info = ensembl.genes_by_name(gene.Name)[0]
        except ValueError:
            name_and_db_not_found.append(gene)
        else:
            gene_names.append(gene.Name)
            db_id.append(gene.Database_ID)
            chromosome.append(info.contig)
            start_bp.append(info.start)
            end_bp.append(info.end)

    # Put what we have so far into a table
    gene_locations = pd.DataFrame({
        "Name": gene_names,
        "Database_ID": db_id,
        "chromosome": chromosome,
        "start_bp": start_bp,
        "end_bp": end_bp
    })

    # Add arms
    cytoband = get_cytoband_info()

    p_arm_max = cytoband[cytoband["arm"] == "p"].\
    groupby("chromosome").\
    max("bp_stop").\
    reset_index(drop=False)[["chromosome", "bp_stop"]].\
    rename(columns={"bp_stop": "p_arm_max"})

    gene_locations = gene_locations.merge(
        p_arm_max,
        on="chromosome",
        how="left",
    )

    gene_locations = gene_locations.assign(
        arm=np.where(gene_locations["start_bp"] <= gene_locations["p_arm_max"], "p", "q")
    )

    # Drop the p_arm_max column now that we have arms
    gene_locations = gene_locations.drop(columns="p_arm_max")

    # Set index as gene name and database ID
    gene_locations = gene_locations.set_index(["Name", "Database_ID"], drop=True, append=False)

    # Package these extra lists
    not_found_genes = {
        "name_and_db_not_found": name_and_db_not_found,
        "release_too_new": release_too_new,
        "old_db_not_found": old_db_not_found,
        "old_db_release_not_found": old_db_release_not_found,
        "no_db_id": no_db_id,
    }

    return gene_locations, not_found_genes

def load_ensembl_release(release_number):
    ensembl = pyensembl.EnsemblRelease(release_number)
    try:
        ensembl.genes() # If this fails, we need to download the data again.
    except ValueError as e:
        print(f"Downloading Ensembl release {release_number} data...", end="\r")
        ensembl.download()
        print(f"Indexing Ensembl release {release_number} data...   ", end="\r")
        ensembl.index()
        print(" " * 80, end="\r") # Clear the message
    return ensembl

def get_cytoband_info():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cytoband_file = os.path.join(BASE_DIR, "cnvutils", "data", 'NCBI_ideogram.csv')
    df = pd.\
    read_csv(cytoband_file).\
    rename(columns={"#chromosome": "chromosome"})
    return df

# Helper functions

def _lookup_old_release(gene_db_id):

    url = f"https://rest.ensembl.org/archive/id/{gene_db_id}"
    params = {"content-type": "application/json"}

    response = requests.get(url, params=params)
    response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful

    return int(response.json()["release"])

def _load_cptac_tables(cancer_types, data_types, pancan, no_internet=False):
    """Get the tables for the specified data types from the specified cancer types.

    Parameters:
    cancer_types (list of str): The cancer types to get data from
    data_types (list of str): The data types to get from each cancer, e.g. proteomics, CNV, transcriptomics, etc.
    pancan (bool): If False, use the regular cptac datasets. If true, use the cptac.pancan (harmonized) datasets.
    no_internet (bool): If True, don't try to update indices. Default False.

    Returns:
    dict of str: dict of str: pd.DataFrame: A dict where the keys are data types and the values are dicts where the keys are cancer types and the values are dataframes of the proper data type.
    """
    # Initialize the dict we'll use to return tables
    all_tables = {}
    for data_type in data_types:
        all_tables[data_type] = {}

    # Load and save tables
    for cancer_type in cancer_types:
        cancer_type_tables = _load_cancer_type_cptac_tables(cancer_type, data_types, pancan, no_internet)
        for data_type, df in cancer_type_tables.items():
            all_tables[data_type][cancer_type] = df

    return all_tables


def _load_cancer_type_cptac_tables(cancer_type, data_types, pancan, no_internet=False):
    """Load the specified data tables from the given cancer type. We have this as a separate function instead of as part of _load_cptac_tables so that the cancer dataset object will be allowed to be garbage collected after we're done with it, instead of sticking around and wasting RAM.

    Parameters:
    cancer_type (str): The cancer type to load
    data_types (list of str): The tables to get
    pancan (bool): If False, use the regular cptac datasets. If true, use the cptac.pancan (harmonized) datasets.
    no_internet (bool): If True, don't try to update indices. Default False.

    Returns:
    dict of str: pd.DataFrame: The requested tables from the given cancer type, indexed by name.
    """

    # Load the cancer type
    if pancan:
        if cancer_type == "brca":
            ds = cptac.pancan.PancanBrca(no_internet=no_internet)
        elif cancer_type == "ccrcc":
            ds = cptac.pancan.PancanCcrcc(no_internet=no_internet)
        elif cancer_type == "coad":
            ds = cptac.pancan.PancanCoad(no_internet=no_internet)
        elif cancer_type == "gbm":
            ds = cptac.pancan.PancanGbm(no_internet=no_internet)
        elif cancer_type == "hnscc":
            ds = cptac.pancan.PancanHnscc(no_internet=no_internet)
        elif cancer_type == "lscc":
            ds = cptac.pancan.PancanLscc(no_internet=no_internet)
        elif cancer_type == "luad":
            ds = cptac.pancan.PancanLuad(no_internet=no_internet)
        elif cancer_type == "ov":
            ds = cptac.pancan.PancanOv(no_internet=no_internet)
        elif cancer_type == "pdac":
            ds = cptac.pancan.PancanPdac(no_internet=no_internet)
        elif cancer_type == "ucec":
            ds = cptac.pancan.PancanUcec(no_internet=no_internet)
        else:
            raise ValueError(f"Invalid cancer type name '{cancer_type}'")

    else:
        if cancer_type == "brca":
            ds = cptac.Brca(no_internet=no_internet)
        elif cancer_type == "ccrcc":
            ds = cptac.Ccrcc(no_internet=no_internet)
        elif cancer_type == "colon":
            ds = cptac.Colon(no_internet=no_internet)
        elif cancer_type == "endometrial":
            ds = cptac.Endometrial(no_internet=no_internet)
        elif cancer_type == "gbm":
            ds = cptac.Gbm(no_internet=no_internet)
        elif cancer_type == "hnscc":
            ds = cptac.Hnscc(no_internet=no_internet)
        elif cancer_type == "lscc":
            ds = cptac.Lscc(no_internet=no_internet)
        elif cancer_type == "luad":
            ds = cptac.Luad(no_internet=no_internet)
        elif cancer_type == "ovarian":
            ds = cptac.Ovarian(no_internet=no_internet)
        elif cancer_type == "pdac":
            ds = cptac.Pdac(no_internet=no_internet)
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
                tables[data_type] = ds.get_transcriptomics(source="washu")
            else:
                raise ValueError(f"Invalid data type name '{data_type}'")
                
        else:
            tables[data_type] = ds._get_dataframe(data_type, tissue_type="both")

    return tables

def _load_gistic_tables(base_dir=os.getcwd())

    gistic_dir = os.path.join(base_dir, "data", "sources", "Broad_pipeline_wxs")
    data_files = [
        "all_lesions.txt",
        "all_data_by_genes.txt",
    ]

    mapping_file_path = os.path.join(base_dir, "data", "sources", "GISTIC_Matched_Samples_Updated.txt")

# Old

def get_driver_genes():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file = os.path.join(BASE_DIR, 'cnvutils', "data", 'bailey_driver_genes.csv')
    df = pd.read_csv(file, skiprows=3)
    return df


def get_normal_expr_table():
    """Load the table of normal protein expression levels for different tissues. This table was downloaded from:
    https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20188503&file=msb188503-sup-0007-TableEV5.zip

    It was produced as part of this paper:
    Wang D, Eraslan B, Wieland T, et al. A deep proteome and transcriptome abundance atlas of 29 healthy human
    tissues. Mol Syst Biol. 2019;15(2):e8503. Published 2019 Feb 18. doi:10.15252/msb.20188503
    """
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(BASE_DIR, "cnvutils", "data", "Table_EV5.xlsx")

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
