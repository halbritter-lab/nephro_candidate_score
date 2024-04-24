# import basic modules
import gzip
import numpy as np
import pandas as pd
import os
import re
import sys
import yaml

# import third-party modules
import vcfpy

# import helper functions
from vs_helper_functions import *

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

# get current date
current_date = datetime.today().strftime('%Y-%m-%d') # Note: problem may arise if Snakemake rule execution takes longer than the same day it was started

# get chromosome
chrom = str(sys.argv[2])
print(f"chromosome: {chrom}")

# get positive genes with evidence count 2-5 from kidney-genetics
pos_genes_dest_file = f"raw/A_MergeAnalysesSources.{config_vars['kidney_genetics_version_vs']}.csv.gz"
pos_genes = pd.read_csv(pos_genes_dest_file, compression='gzip')
pos_genes2345 = pos_genes.query("evidence_count > 1")
pos_genes2345_hgnc_id_list = pos_genes2345['hgnc_id'].tolist()


# get latest VEP annotated VCF
latest_file = get_latest_file(prefix="clinvar_vars_kid-gen_2345_vep_anno", directory="features_labels/results/", extension=None) 
print(f"Extract variant features from '{latest_file}'.")

# open VEP annotated VCF, this will read in the header
reader = vcfpy.Reader.from_path(latest_file)

# check if VCF has a FILTER line in header
if len(reader.header.get_lines("FILTER").mapping.keys()) > 0:
    raise ValueError("FILTERS are present. Currently not accounted for in code.") # TODO: change??

# get VEP 'CSQ' annotation 
CSQ_info = reader.header.get_info_field_info('CSQ')
CSQ_description = CSQ_info.description
CSQ_description

# extract pattern word|word|...|word for getting the field names
pattern = r'\b(?:\w+\|)+\w+\b'
matches = re.findall(pattern, CSQ_description)

CSQ_fields = matches[0].split('|')


# choose CSQ fields
selected_fields = ['Consequence', 'IMPACT', 'HGNC_ID', 'SYMBOL', 'Feature', 'gnomADe_AF', 'gnomADg_AF', 'CADD_PHRED', 'PUBMED', 'NMD', 'BLOSUM62', 'MaxEntScan_alt', 'MaxEntScan_diff', 'MaxEntScan_ref', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'PrimateAI_rankscore', 'REVEL_rankscore', 'phastCons100way_vertebrate_rankscore']

# create an empty df
clinvar_vep_anno = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'type', 'QUAL', 'CLNSIG'] + selected_fields)

count = 0 

# read VCF
for record in reader:
    if str(record.CHROM) == chrom:

        # get basic info of each record: CHROM, POS, REF, ALT, ...
        base_df = pd.DataFrame({
            'CHROM' : record.CHROM,
            'POS' : record.POS,
            'ID' : record.ID,
            'REF': record.REF,
            'ALT': [i.value for i in record.ALT],
            'type': [i.type for i in record.ALT],
            'QUAL': record.QUAL,
            'CLNSIG': [i for i in record.INFO['CLNSIG']],
            'dummy_key': 1 # needed for join later
        })

        # get VEP annotation values 'CSQ'
        INFO = record.INFO
        CSQ_cell_entries = [i.split("|") for i in INFO['CSQ']]
        CSQ_record_df = pd.DataFrame(CSQ_cell_entries, columns=CSQ_fields)
        CSQ_record_df['dummy_key'] = 1 # needed for join

        # merge CSQ df with basic info df
        new_rows = pd.merge(base_df, CSQ_record_df, on='dummy_key')[list(base_df.columns) + selected_fields]

        # concatenate summarizing df with new entries 
        clinvar_vep_anno = pd.concat([clinvar_vep_anno, new_rows], axis=0)

    # increase count to monitor progress
    count = count + 1
    if count % 10000 == 0: 
        print(count)
        sys.stdout.flush() 

        
print("All records read.")
sys.stdout.flush() 

clinvar_vep_anno['hgnc_id_int'] = clinvar_vep_anno['HGNC_ID'].str.split(":").str[1]

# filter out NaN values for HGNC ID
clinvar_vep_anno = clinvar_vep_anno.query("hgnc_id_int.notna()")

# cast HGNC id to integer
clinvar_vep_anno['hgnc_id_int'] = pd.to_numeric(clinvar_vep_anno['hgnc_id_int'], downcast='integer')

# filter out genes that are not in kidney-genetics (evidence count 2-5)
clinvar_kid_gen2345_vep_anno = clinvar_vep_anno.query('hgnc_id_int in @pos_genes2345_hgnc_id_list')

# write csv
clinvar_kid_gen2345_vep_anno.to_csv(f"features_labels/results/raw_features_clinvar_vars_kid-gen_2345_chr{chrom}_{current_date}.csv.gz", index=False, compression='gzip') 
