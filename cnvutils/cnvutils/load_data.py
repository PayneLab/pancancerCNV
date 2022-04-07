import collections
import cptac
import cptac.pancan
import glob
import json
import numpy as np
import os
import pandas as pd
import pyensembl
import re
import requests
import sys
import xmltodict

from .constants import ALL_CANCERS

def save_input_tables(pancan, data_dir=os.path.join(os.getcwd(), "..", "data")):
    """Load CNV, transcriptomics, and proteomics tables for all cancers, create
    gene locations table, and save all of them.
    """

    # Load cptac tables
    cptac_tables = _load_cptac_tables(
        cancer_types=ALL_CANCERS[0],
        data_types=[
            "CNV",
            "transcriptomics",
            "proteomics",
        ],
        pancan=pancan,
    )

    # Load, reformat, and save the GISTIC tables
    gistic_tables = _load_gistic_tables()

    # Get list of all genes we have CNV data for. These are the ones we'll need the locations of.
    genes = None
    for df in cptac_tables["CNV"].values() + gistic_tables["gene"].values():
        if genes is None:
            genes = df.columns.copy(deep=True) # Just for the first one
        else:
            genes = genes.union(df.columns)

    genes = genes.to_frame()
    gene_locations, not_found_genes = query_gene_locations_database(genes)

    # TODO
    # Do some comparisons of the Entrez and Ensembl gene location data

    # Create a data directory in the directory the function was called from
    input_data_dir = os.path.join(data_dir, "sources")
    cptac_data_dir = os.path.join(input_data_dir, "cptac_tables")
    gistic_data_dir = os.path.join(input_data_dir, "gistic_tables_formatted")

    os.makedirs(input_data_dir, exist_ok=True)
    os.makedirs(cptac_data_dir, exist_ok=True)
    os.makedirs(gistic_data_dir, exist_ok=True)

    # Save the gene locations
    gene_locations_save_path = os.path.join(input_data_dir, "gene_locations.tsv.gz")
    gene_locations.to_csv(gene_locations_save_path, sep="\t")

    # Save the cptac omics tables
    for data_type, cancer_types in cptac_tables.items():
        for cancer_type, df in cancer_types.items():

            file_name = f"{cancer_type}_{data_type}.tsv.gz"
            save_path = os.path.join(cptac_data_dir, file_name)
            df.to_csv(save_path, sep="\t")

    # Save GISTIC CNV tables
    for data_level, cancer_types in gistic_tables.items():
        for cancer_type, df in cancer_types.items():

            file_name = f"{cancer_type}_{data_level}.tsv.gz"
            save_path = os.path.join(gistic_data_dir, file_name)
            df.to_csv(save_path, sep="\t")

def load_input_tables(data_dir, data_types=["CNV", "proteomics", "transcriptomics"], cancer_types=ALL_CANCERS[0]):

    # Standardize data_types and cancer_types
    data_types = [data_type.lower() for data_type in data_types]
    cancer_types = [cancer_type.lower() for cancer_type in cancer_types]

    # Get the data tables directory
    cptac_tables_dir = os.path.join(data_dir, "sources", "cptac_tables")

    # Get list of tables to load
    all_table_paths = sorted(glob.glob(os.path.join(cptac_tables_dir, "*")))
    load_table_paths = []
    for path in all_table_paths:
        cancer_type, data_type = path.split(os.sep)[-1].split(".")[0].split("_")
        if cancer_type.lower() in cancer_types and data_type.lower() in data_types:
            load_table_paths.append(path)

    # Load the tables
    cptac_tables = {}
    for i, path in enumerate(load_table_paths):
        cancer_type, data_type = path.split(os.sep)[-1].split(".")[0].split("_")

        print(f"Loading {cancer_type} {data_type} ({i + 1}/{len(load_table_paths)})...{' ' * 30}", end="\r")

        if data_type not in cptac_tables.keys():
            cptac_tables[data_type] = {}

        cptac_tables[data_type][cancer_type] = pd.read_csv(
            path,
            sep="\t",
            index_col=0,
            header=[0, 1]
        )

    # Clear last loading message
    print(" " * 80, end="\r")

    return cptac_tables

def load_gene_locations(data_dir=os.path.join(os.getcwd(), "..", "data")):

    gene_locations_path = os.path.join(data_dir, "sources", "gene_locations.tsv.gz")
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

def _lookup_old_ensembl_release(gene_db_id):

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

def _load_gistic_tables(levels, data_dir=os.path.join(os.getcwd(), "..", "data")):

    # Load and format mapping file
    mapping_file_path = os.path.join(data_dir, "sources", "GISTIC_Matched_Samples_Updated.txt")

    id_map = pd.\
    read_csv(mapping_file_path, sep="\t").\
    rename(columns={"Case_ID": "Patient_ID"}).\
    melt(
        id_vars=["Patient_ID", "Tumor_Type", "Sample_Type", "Category"],
        var_name="other_id_type",
        value_name="other_id",
        ignore_index=True,
    ).\
    drop(columns=["Sample_Type", "Category", "other_id_type"]).\
    drop_duplicates(keep="first")

    # Add an empty level to the columns index to match the format of the data table
    id_map = id_map.\
    transpose().\
    assign(NCBI_ID=[""] * id_map.shape[1]).\
    set_index("NCBI_ID", append=True).\
    transpose()

    # Set the column multiindex level names
    id_map.columns.names = ["Name", "NCBI_ID"]

    # Get the requested data types
    gistic_dir = os.path.join(data_dir, "sources", "Broad_pipeline_wxs")
    data_file_paths = {
        "segment_level": os.path.join(gistic_dir, "all_lesions.txt"),
        "gene_level": os.path.join(gistic_dir, "all_data_by_genes.txt"),
        "arm_level": os.path.join(gistic_dir, "broad_values_by_arm.txt"),
    }

    tables = {}
    for level in levels:

        if level == "segment":
            df = pd.read_csv(data_file_paths["segment_level"], sep="\t")

            # Drop empty last row. It must have been a formatting error when the file was generated.
            df = df.drop(columns="Unnamed: 1092")

            # Select the rows with raw CNV values, not thresholded
            df = df[df["Unique Name"].str.endswith("values")]

            df = df.drop(columns=[
                "Unique Name",
                "Descriptor",
                "q values",
                "Residual q values after removing segments shared with higher peaks",
                "Broad or Focal",
                "Amplitude Threshold",
                #"Wide Peak Limits", # TODO: Which limits to use? Region?
                #"Peak Limits",
                #"Region Limits",
            ])

            import pdb; pdb.set_trace()
            print(1)

        elif level == "gene":
            df = pd.read_csv(data_file_paths["gene_level"], sep="\t").\
            drop(columns="Cytoband").\
            rename(columns={
                "Gene Symbol": "Name",
                "Gene ID": "NCBI_ID",
            })

            # Some NCBI IDs are negative, which are incorrect. Query Entrez with the locus name to get the correct ID, then put those in.
            neg_ids = df[df["NCBI_ID"] < 0]
            neg_ids = neg_ids.assign(correct_id=neg_ids["Name"].apply(_lookup_ncbi_id_by_name))
            df.loc[neg_ids.index, "NCBI_ID"] = neg_ids["correct_id"]

            # Get gene metadata from NCBI Entrez Gene database
            gene_ids = df["NCBI_ID"].drop_duplicates(keep="first").astype(str)
            metadata = _lookup_genes_ncbi_id(gene_ids)

            # Look up data for new IDs for any genes listed as "secondary" (meaning they've been replaced by another gene)
            secondary = metadata[metadata["status"] == "secondary"]
            new_meta = _lookup_genes_ncbi_id(secondary["new_id"])

            # TODO
            # Do we drop discontinued genes that still have start and end bp?
            # For genes with the wrong names, replace the names in df with the new names
            # For secondary genes, replace the IDs ad names in df with the new names and IDs from new_meta
            # Split into a CNV data table and a location table

            import pdb; pdb.set_trace()

            # (Pdb) live.shape
            # (25739, 8)

            # (Pdb) sec.shape
            # (74, 8)

            # (Pdb) dead.shape
            # (45, 8)

            # (Pdb) new_meta.isna().sum()
            # Name           0
            # chromosome     0
            # status         0
            # new_id        74
            # NCBI_ID        0
            # Ensembl_ID    11
            # start_bp       0
            # end_bp         0
            # dtype: int64

            # (Pdb) metadata.isna().sum()
            # Name              0
            # chromosome        0
            # status            0
            # new_id        25784
            # NCBI_ID           0
            # Ensembl_ID     1896
            # start_bp         95
            # end_bp           95
            # dtype: int64
            
            # (Pdb) live.isna().sum()
            # Name              0
            # chromosome        0
            # status            0
            # new_id        25739
            # NCBI_ID           0
            # Ensembl_ID     1777
            # start_bp         16
            # end_bp           16
            # dtype: int64

            # (Pdb) sec.isna().sum()
            # Name           0
            # chromosome     0
            # status         0
            # new_id         0
            # NCBI_ID        0
            # Ensembl_ID    74
            # start_bp      72
            # end_bp        72
            # dtype: int64
            
            # (Pdb) dead.isna().sum()
            # Name           0
            # chromosome     0
            # status         0
            # new_id        45
            # NCBI_ID        0
            # Ensembl_ID    45
            # start_bp       7
            # end_bp         7
            # dtype: int64

            # Resolve duplicate versions of different genes shown on different chromosomes

            # Join in Ensembl identifiers

            # TODO: When the Broad gives us the new tables without duplicated genes or negative IDs, see if we still have any duplicated NCBI IDs. Hopefully we won't.

            # Save the gene metadata for later
            tables["gene_metadata"] = metadata
            
            # Set index columns and transpose
            df = df.set_index(["Name", "NCBI_ID"]).\
            transpose()

        elif level == "arm":
            pass

        else:
            raise ValueError(f"Invalid GISTIC data level: '{level}'")

        # Merge in mapping index
        df = df.\
        merge(
            right=id_map,
            how="left",
            left_index=True,
            right_on="other_id",
        ).\
        drop(columns=("other_id", ""))

        # Standardize cancer types
        df = df.\
        assign(Tumor_Type=df["Tumor_Type"].replace({
            "BR": "brca",
            "CCRCC": "ccrcc",
            "CO": "coad",
            "GBM": "gbm",
            "HNSCC": "hnscc",
            "LSCC": "lscc",
            "LUAD": "luad",
            "OV": "ov",
            "PDA": "pdac",
            "UCEC": "ucec",
        })).\
        rename(columns={"Tumor_Type": "cancer_type"}).\
        sort_values(by=["cancer_type", "Patient_ID"]).\
        set_index(["Patient_ID", "cancer_type"])

        # Split into a separate table for each cancer type
        cancer_types_dict = {cancer_type: df for cancer_type, df in df.groupby("cancer_type")}

        tables[level] = cancer_types_dict

    return tables

def _lookup_ncbi_id_by_name(name, field="title"):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    # Run the query
    params = {
        "db": "gene",
        "term": name,
        "field": field,
        "retmode": "xml",
        "api_key": _get_entrez_api_key()
    }

    response = requests.get(url, params=params)
    response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful

    # Parse XML string as a nested dictionary
    resp_xml = xmltodict.parse(response.text)
    if resp_xml["eSearchResult"]["Count"] == "0":
        if field != "all":
            ncbi_id = _lookup_ncbi_id_by_name(name, field="all")
        else:
            ncbi_id = np.nan
    elif isinstance(resp_xml["eSearchResult"]["IdList"]["Id"], list):
        ncbi_id = int(resp_xml["eSearchResult"]["IdList"]["Id"][0])
    else:
        ncbi_id = int(resp_xml["eSearchResult"]["IdList"]["Id"])

    # Return the ID
    return ncbi_id

def _lookup_genes_ncbi_id(gene_ids):

    # The NCBI Entrez API won't except query strings that are too long. For this particular query,
    # the maximum character length of the gene ID portion of the query string is 3122 characters.
    # So, we'll split our query into batches based on that length.
    escaped_comma = "%2C" # This is the URL escape character for a comma
    full_ids_str = escaped_comma.join(gene_ids)
    max_chars = 4040 # Specifically, this is the maximum number of characters possible for the id parameter of the query URL, while still leaving enough room for the rest of the URL. So the total length of the URL will be a little longer than this. 

    queried = []

    all_genes_info = pd.DataFrame()
    while len(full_ids_str) > 0:

        # See if what's left is too long
        if len(full_ids_str) > max_chars:

            # Find the end of the last ID before the query string would be too long
            last_comma_index = max_chars
            while full_ids_str[last_comma_index:last_comma_index + len(escaped_comma)] != escaped_comma:
                last_comma_index -= 1

            # Slice out that much of the query string, and leave the rest for next time, minus the leading comma
            ids_str_slice = full_ids_str[:last_comma_index]
            full_ids_str = full_ids_str[last_comma_index + len(escaped_comma):]

        # If it's not to long, just use the whole thing
        else:
            ids_str_slice = full_ids_str
            full_ids_str = ""

        # Print an info message about what we're looking up
        total_genes = len(gene_ids)
        end = total_genes - len(full_ids_str.split(escaped_comma))
        start = end - len(ids_str_slice.split(escaped_comma))
        print(f"Looking up genes {start + 1} to {end} of {total_genes}...", end="\r")

        # ids_str_slice currently uses %2C as the separator value between ids, which is the URL encoded hexadecimal value 
        # for a comma. However, the requests library automatically tries to encode the URL for us, and won't recognize this
        # as an already escaped character. So, it will see the "%" as a standalone character and try to escape that as %25,
        # which will then mess up our escaped commas by making them read %252C. So, to avoid this double encoding problem,
        # we're going to replace the %2C values with regular commas before we send the string through the requests library.
        # You might ask why we didn't just use regular commas in the first place, and not worry about escaping them manually
        # if that won't work anyways. The reason is because the Entrez API limits the length of URLs it will accept requests
        # through, and to test the actual length of the URL so we can know how many gene IDs we can send at once, we need to
        # use the escape sequence that regular commas will be replaced with. So up to this point we've used the %2C escape
        # sequence, but now we'll replace that with commas before sending the string to the requests.get call.
        comma = ","
        ids_str_slice = comma.join(ids_str_slice.split(escaped_comma))

        # Run the query in a separate function so that any variables not returned (including all the text we don't need from
        # the response body) will be garbage collected when the function ends. This will save RAM.
        genes_info = _run_ncbi_id_query(ids_str_slice)

        # TEMP: check
        queried.extend(ids_str_slice.split(","))

        # Save the results
        all_genes_info = pd.concat(
            [all_genes_info, genes_info],
            axis="index"
        )

        # Clear the info message
        print(" " * 100, end="\r")

    return all_genes_info

def _run_ncbi_id_query(ids_str_slice):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    # Run the query
    params = {
        "db": "gene",
        "id": ids_str_slice,
        "retmode": "xml",
    }

    response = requests.get(url, params=params)
    response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful

    # Parse XML string as a nested dictionary
    resp_xml = xmltodict.parse(response.text)

    # Extract needed info for each gene
    genes = []
    ncbi_ids = []
    statuses = []
    new_ids = []
    ensembl_ids = []
    chrs = []
    starts_bp = []
    ends_bp = []
    for gene_xml in resp_xml["Entrezgene-Set"]["Entrezgene"]:

        # Get the gene
        gene = gene_xml["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]

        # Get the NCBI Entrez gene ID
        ncbi_id = gene_xml["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]

        # Get the identifier status
        status = gene_xml["Entrezgene_track-info"]["Gene-track"]["Gene-track_status"]["@value"] 

        # See if there's a different current ID
        new_id = np.nan 
        if "Gene-track_current-id" in gene_xml["Entrezgene_track-info"]["Gene-track"].keys():
            if isinstance(gene_xml["Entrezgene_track-info"]["Gene-track"]["Gene-track_current-id"]["Dbtag"], list):
                for db in gene_xml["Entrezgene_track-info"]["Gene-track"]["Gene-track_current-id"]["Dbtag"]:
                    if db["Dbtag_db"] == "GeneID":
                        new_id = db["Dbtag_tag"]["Object-id"]["Object-id_id"]
                        break

        # Get the Ensembl ID
        ensembl_id = np.nan
        if "Gene-ref_db" in gene_xml["Entrezgene_gene"]["Gene-ref"].keys():
            if isinstance(gene_xml["Entrezgene_gene"]["Gene-ref"]["Gene-ref_db"]["Dbtag"], list):
                for db_info in gene_xml["Entrezgene_gene"]["Gene-ref"]["Gene-ref_db"]["Dbtag"]:
                    if db_info["Dbtag_db"] == "Ensembl":
                        ensembl_id = db_info["Dbtag_tag"]["Object-id"]["Object-id_str"]
                        break
            else:
                db_info = gene_xml["Entrezgene_gene"]["Gene-ref"]["Gene-ref_db"]["Dbtag"]
                if db_info["Dbtag_db"] == "Ensembl":
                    ensembl_id = db_info["Dbtag_tag"]["Object-id"]["Object-id_str"]

        # Get the chromosome
        if "BioSource_subtype" in gene_xml["Entrezgene_source"]["BioSource"].keys():
            chrm = gene_xml["Entrezgene_source"]["BioSource"]["BioSource_subtype"]["SubSource"]["SubSource_name"]
        else:
            chrm = np.nan

        # Get the start and end base pairs
        start_bp = np.nan
        end_bp = np.nan
        if isinstance(gene_xml["Entrezgene_locus"]["Gene-commentary"], list):
            for comment in gene_xml["Entrezgene_locus"]["Gene-commentary"]:
                if "Gene-commentary_heading" in comment.keys() and comment["Gene-commentary_heading"].endswith("Primary Assembly"):
                    start_bp = comment["Gene-commentary_seqs"]["Seq-loc"]["Seq-loc_int"]["Seq-interval"]["Seq-interval_from"]
                    end_bp = comment["Gene-commentary_seqs"]["Seq-loc"]["Seq-loc_int"]["Seq-interval"]["Seq-interval_to"]
                    break
        elif "Gene-commentary_seqs" in gene_xml["Entrezgene_locus"]["Gene-commentary"].keys():
            start_bp = gene_xml["Entrezgene_locus"]["Gene-commentary"]["Gene-commentary_seqs"]["Seq-loc"]["Seq-loc_int"]["Seq-interval"]["Seq-interval_from"]
            end_bp = gene_xml["Entrezgene_locus"]["Gene-commentary"]["Gene-commentary_seqs"]["Seq-loc"]["Seq-loc_int"]["Seq-interval"]["Seq-interval_to"]

        # Save the info
        genes.append(gene)
        ncbi_ids.append(ncbi_id)
        statuses.append(status)
        new_ids.append(new_id)
        ensembl_ids.append(ensembl_id)
        chrs.append(chrm)
        starts_bp.append(start_bp)
        ends_bp.append(end_bp)

    results = pd.DataFrame({
        "Name": genes,
        "chromosome": chrs,
        "status": statuses,
        "new_id": new_ids,
        "NCBI_ID": ncbi_ids,
        "Ensembl_ID": ensembl_ids,
        "start_bp": starts_bp,
        "end_bp": ends_bp,
    })

    return results

def _get_entrez_api_key():

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    api_key_path = os.path.join(base_dir, "cnvutils", "entrez_api_key.txt")

    if not os.path.isfile(api_key_path):
        raise ValueError(f"You must create an Entrez API key. You can do this by going to <https://www.ncbi.nlm.nih.gov/account/>, creating an account, and then navigating to NCBI Site Preferences>Account Settings. Then, take the API key and save it in a plaintext file located at '{api_key_path}'.")

    with open(api_key_path, "r") as file:
        api_key = file.read().strip()

    return api_key

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
