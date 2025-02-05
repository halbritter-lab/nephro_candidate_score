# import basic modules
import gzip
import numpy as np
import pandas as pd
import os
import urllib.request
import yaml

# import third-party modules
import vcfpy

# get global config file
CONFIG_FILE = os.getenv('CONFIG_FILE')

# define relative script path
project_topic = "nephrology"
project_name = "nephro_candidate_score"
script_path = "/variant_score/"

# read configs
with open(CONFIG_FILE, 'r') as file:
    config_data = yaml.safe_load(file)

config_vars = config_data[project_topic]

# set working directory
os.chdir(f"{config_vars['ML_projectsdir']}{project_name}{script_path}")


# download HGNC annotated table from kidney-genetics
hgnc_annotated_url = f"https://raw.githubusercontent.com/halbritter-lab/kidney-genetics/main/analyses/B_AnnotationHGNC/results/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version_vs']}.csv.gz"
hgnc_annotated_dest_file = f"raw/non_alt_loci_set_coordinates.{config_vars['hgnc_gt_version_vs']}.csv.gz"

# check if the file already exists
if not os.path.exists(hgnc_annotated_dest_file):
    # download the file
    urllib.request.urlretrieve(hgnc_annotated_url, hgnc_annotated_dest_file)
    print(f"The file '{hgnc_annotated_dest_file}' has been downloaded.")
else:
    print(f"The file '{hgnc_annotated_dest_file}' already exists. Skipping the download.")

# read in file    
hgnc_annotated = pd.read_csv(hgnc_annotated_dest_file, compression='gzip', low_memory=False)

# add a new column without the "HGNC:" prefix
hgnc_annotated['hgnc_id_int'] = hgnc_annotated['hgnc_id'].str.replace('HGNC:', '')

# convert the 'hgnc_id_int' and 'entrez_id' column to integers
hgnc_annotated['hgnc_id_int'] = pd.to_numeric(hgnc_annotated['hgnc_id_int'], downcast='integer')

# download positive genes from kidney-genetics
pos_genes_url = f"https://github.com/halbritter-lab/kidney-genetics/raw/main/analyses/A_MergeAnalysesSources/results/A_MergeAnalysesSources.{config_vars['kidney_genetics_version_vs']}.csv.gz"
pos_genes_dest_file = f"raw/A_MergeAnalysesSources.{config_vars['kidney_genetics_version_vs']}.csv.gz"

# check if the file already exists
if not os.path.exists(pos_genes_dest_file):
    # download the file
    urllib.request.urlretrieve(pos_genes_url, pos_genes_dest_file)
    print(f"The file  '{pos_genes_dest_file}' has been downloaded.")
else:
    print(f"The file '{pos_genes_dest_file}' already exists. Skipping the download.")

# read in file    
pos_genes = pd.read_csv(pos_genes_dest_file, compression='gzip')

# filter only genes with evidence count 2-5, merge with entrez ID from HGNC table, convert entrez ID to integer
pos_genes2345 = pos_genes.query("evidence_count > 1").merge(hgnc_annotated[['hgnc_id_int', 'entrez_id']], how='left', left_on='hgnc_id', right_on='hgnc_id_int')
pos_genes2345['entrez_id'] = pd.to_numeric(pos_genes2345['entrez_id'], downcast='integer')
pos_genes2345_entrez_id_list = pos_genes2345['entrez_id'].tolist()



## Filter Clinvar VCF file for variants in kidney-genetics genes with evidence count 2-5
# open VCF, this will read in the header
reader = vcfpy.Reader.from_path(f"raw/clinvar_{config_vars['clinvar_version']}.vcf.gz")

# create a VCF writer
writer = vcfpy.Writer.from_path('/dev/stdout', reader.header)

# read VCF
for record in reader:
    # check if GENEINFO is available per variant    
    if 'GENEINFO' in record.INFO.keys():
        GENEINFO = record.INFO['GENEINFO']
        GENEINFO_entrez_ids = [int(i.split(":")[1]) for i in GENEINFO.split("|")]

    # if variant is in one of the kidney-genetics genes write it to subset VCF
        if any(elem in pos_genes2345_entrez_id_list for elem in GENEINFO_entrez_ids):
            writer.write_record(record)

            
# close writer    
writer.close()