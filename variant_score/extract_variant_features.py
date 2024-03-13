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

# get positive genes with evidence count 2-5 from kidney-genetics
pos_genes_dest_file = f"raw/A_MergeAnalysesSources.{config_vars['kidney_genetics_version_vs']}.csv.gz"
pos_genes = pd.read_csv(pos_genes_dest_file, compression='gzip')
pos_genes2345 = pos_genes.query("evidence_count > 1")
pos_genes2345_hgnc_id_list = pos_genes2345['hgnc_id'].tolist()




# get latest VEP annotated VCF
latest_file = get_latest_file(prefix="clinvar_vars_kid-gen_2345_vep_anno", directory=".", extension=None) #TODO: change!
print(f"Extract variant features from '{latest_file}'.")

# open VEP annotated VCF, this will read in the header
# reader = vcfpy.Reader.from_path('clinvar_vars_kid-gen_2345_vep_anno_2024-03-11.vcf.gz')
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
CSQ_df = pd.DataFrame(columns=CSQ_fields)

# selected VEP fields
# selected_fields = ['Consequence', 'IMPACT', 'HGNC_ID', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 
#                   'BIOTYPE', 'EXON', 'SIFT', 'PolyPhen', 'gnomADe_AF', 'gnomADg_AF', 'CLIN_SIG', 
#                   'CADD_PHRED', 'CADD_RAW']

# selected_fields = ['Consequence', 'IMPACT', 'HGNC_ID', 'SYMBOL', 'gnomADe_AF', 'gnomADg_AF', 
#                   'CADD_PHRED', 'CADD_RAW']#, #'PUBMED']#, 
#                 #  'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']

    
# selected_fields = ['Consequence', 'IMPACT', 'HGNC_ID', 'SYMBOL', 'gnomADe_AF', 'gnomADg_AF', 
#                    'CADD_PHRED', 'CADD_RAW', 'PUBMED', 
#                    'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']    
    
selected_fields = ['Consequence', 'IMPACT', 'HGNC_ID', 'SYMBOL', 'Feature',
                   'gnomADe_AF', 'gnomADg_AF', 
                  'CADD_PHRED', 
#                    'CADD_RAW', 
                   'PUBMED', 
                   'NMD', 
                   'BLOSUM62', 
                   'MaxEntScan_alt', 'MaxEntScan_diff', 'MaxEntScan_ref',
                 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL',
#                   'PrimateAI_pred', 
                   'PrimateAI_rankscore', 
#                    'PrimateAI_score', 
                   'REVEL_rankscore', 
#                    'REVEL_score', 
#                    'phastCons100way_vertebrate', 
                   'phastCons100way_vertebrate_rankscore'
                  ]

clinvar_vep_anno = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'type', 'QUAL', 'CLNSIG'] + selected_fields)

count = 0 #TODO: remove

# read VCF
for record in reader:
    # get basic info: CHROM, POS, REF, ALT, ...
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
    CSQ_record_df['dummy_key'] = 1
    
    # merge CSQ df with basic info df
    new_rows = pd.merge(base_df, CSQ_record_df, on='dummy_key')[list(base_df.columns) + selected_fields]

    # concatenate summarizing df with new entries 
    clinvar_vep_anno = pd.concat([clinvar_vep_anno, new_rows], axis=0)

#     if clinvar_vep_anno.shape[0] > 100000:
#         break # TODO: remove
        
        
    count = count + 1
    if count % 10000 == 0:
        print(count)
        sys.stdout.flush() 


clinvar_vep_anno['hgnc_id_int'] = clinvar_vep_anno['HGNC_ID'].str.split(":").str[1]

# filter out NaN values for HGNC ID
clinvar_vep_anno = clinvar_vep_anno.query("hgnc_id_int.notna()")

# cast HGNC id to integer
clinvar_vep_anno['hgnc_id_int'] = pd.to_numeric(clinvar_vep_anno['hgnc_id_int'], downcast='integer')

# write csv
clinvar_vep_anno.to_csv(f"raw_features_clinvar_vars_kid-gen_2345_prefiltered_{datetime.today().strftime('%Y-%m-%d')}.csv.gz", index=False, compression='gzip')


# filter out genes that are not in kidney-genetics (evidence count 2-5)
clinvar_kid_gen2345_vep_anno = clinvar_vep_anno.query('hgnc_id_int in @pos_genes2345_hgnc_id_list')

# write csv
clinvar_kid_gen2345_vep_anno.to_csv(f"raw_features_clinvar_vars_kid-gen_2345_{datetime.today().strftime('%Y-%m-%d')}.csv.gz", index=False, compression='gzip') 



# # create an own variant ID
# clinvar_kid_gen2345_vep_anno = clinvar_kid_gen2345_vep_anno.copy()
# clinvar_kid_gen2345_vep_anno.loc[:, 'var_ID'] = clinvar_kid_gen2345_vep_anno['CHROM'].astype(str) + '_' \
# + clinvar_kid_gen2345_vep_anno['POS'].astype(str) + '_' + clinvar_kid_gen2345_vep_anno['REF'].astype(str) \
# + '_' + clinvar_kid_gen2345_vep_anno['ALT'].astype(str)

# # map IMPACT to integer
# impact_mapping = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
# clinvar_kid_gen2345_vep_anno['IMPACT_num'] = clinvar_kid_gen2345_vep_anno['IMPACT'].map(impact_mapping)

        


