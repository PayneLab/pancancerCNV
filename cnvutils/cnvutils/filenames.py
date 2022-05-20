import os

# Input files
def get_input_data_dir(data_dir):
    path = os.path.join(data_dir, "sources")
    os.makedirs(path, exist_ok=True)
    return path

def get_input_source_data_dir(data_dir, source):
    path = os.path.join(get_input_data_dir(data_dir), f"{source}_tables")
    os.makedirs(path, exist_ok=True)
    return path

def get_input_source_file_path(data_dir, source, cancer_type, data_type):
    return os.path.join(get_input_source_data_dir(data_dir, source), f"{cancer_type}_{data_type}.tsv.gz")

# Gene locations
def get_gene_locations_path(data_dir, source):
    return os.path.join(get_input_data_dir(data_dir), f"{source}_gene_locations.tsv.gz")

def get_gene_name_updates_path(data_dir, source):
    return os.path.join(get_input_data_dir(data_dir), f"{source}_gene_name_updates.tsv.gz")

# Analysis files
def get_cnv_counts_path(data_dir, source, level, chromosome):
    return os.path.join(data_dir, f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}_cnv_counts.tsv")

def get_has_event_path(data_dir, source, cancer_type, level, chromosome, arm, gain_or_loss):
    return os.path.join(data_dir, f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}{arm}_{gain_or_loss}_{cancer_type}_has_event.tsv")

def get_ttest_results_path(data_dir, source, level, chromosome, arm, gain_or_loss, cis_or_trans, proteomics_or_transcriptomics):
    return os.path.join(data_dir, f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}{arm}_{gain_or_loss}_{cis_or_trans}_{proteomics_or_transcriptomics}_ttest.tsv")

# Metadata
def get_event_metadata_path(data_dir, source, level, chromosome, arm, gain_or_loss):
    return os.path.join(data_dir, f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}{arm}_{gain_or_loss}_metadata.json")

# Plots
def get_charts_img_path(data_dir):
    path = os.path.join(data_dir, "charts_img")
    os.makedirs(path, exist_ok=True)
    return path

def get_chr_line_plot_path(data_dir, source, level, chromosome, chart_format):
    return os.path.join(get_charts_img_path(data_dir), f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}_line_plot.{chart_format}")

def get_chr_gradient_plot_path(data_dir, source, level, chromosome, chart_format):
    return os.path.join(get_charts_img_path(data_dir), f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}_gradient_plot.{chart_format}")

def get_ttest_plot_path(data_dir, source, level, chromosome, arm, gain_or_loss, cis_or_trans, proteomics_or_transcriptomics, chart_format):
    return os.path.join(get_charts_img_path(data_dir), f"{source}_{level + '_' if level else ''}chr{chromosome:0>2}{arm}_{gain_or_loss}_{cis_or_trans}_{proteomics_or_transcriptomics}_ttest_plot.{chart_format}")
