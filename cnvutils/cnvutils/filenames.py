import os

# Input files
def get_input_data_dir(data_dir):
    return os.path.join(data_dir, "sources")

def get_input_source_data_dir(data_dir, source):
    return os.path.join(get_input_data_dir(data_dir), f"{source}_tables")

def get_input_source_file_path(data_dir, source, cancer_type, data_type):
    return os.path.join(get_input_source_data_dir(data_dir, source), f"{cancer_type}_{data_type}.tsv.gz")

# Gene locations
def get_gene_locations_path(data_dir, source):
    return os.path.join(data_dir, "sources", f"{source}_gene_locations.tsv.gz")

def get_gene_name_changes_path(data_dir, source):
    return os.path.join(data_dir, "sources", f"{source}_gene_name_updates.tsv.gz")

# Analysis files
def get_cnv_counts_path(data_dir, source, level, chromosome):
    return os.path.join(data_dir, f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}_cnv_counts.tsv")

def get_has_event_path(data_dir, source, cancer_type, level, chromosome, arm, gain_or_loss):
    return os.path.join(data_dir, f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}{arm}_{gain_or_loss}_{cancer_type}_has_event.tsv")

def get_ttest_results_path(data_dir, source, level, chromosome, arm, gain_or_loss, cis_or_trans):
    return os.path.join(data_dir, f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}{arm}_{gain_or_loss}_{cis_or_trans}_ttest.tsv")
