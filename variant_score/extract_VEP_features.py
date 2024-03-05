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



# open VEP annotated VCF, this will read in the header
reader = vcfpy.Reader.from_path('./clinvar_vars_kid-gen_2345_vep_anno_2024-03-05.vcf.gz')

# get VEP 'CSQ' annotation 
CSQ_info = reader.header.get_info_field_info('CSQ')
CSQ_description = CSQ_info.description

# extract pattern word|word|...|word for getting the 'CSQ' field names
pattern = r'\b(?:\w+\|)+\w+\b'
matches = re.findall(pattern, CSQ_description)
CSQ_fields = matches[0].split('|')


# select VEP 'CSQ' fields
selected_fields = ['Consequence', 'IMPACT', 'HGNC_ID', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 
                  'BIOTYPE', 'EXON', 'SIFT', 'PolyPhen', 'gnomADe_AF', 'gnomADg_AF', 'CLIN_SIG', 
                  'CADD_PHRED', 'CADD_RAW']

# create empty dataframes
CSQ_df = pd.DataFrame(columns=CSQ_fields)
clinvar_kid_gen2345_vep_anno = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'type', 'QUAL', 'FILTER'] + selected_fields)



count = 0

# read VCF
for record in reader:
    
    # get basic info CHROM, POS, REF, ALT, ...
    base_df = pd.DataFrame({
        'CHROM' : record.CHROM,
        'POS' : record.POS,
        'ID' : record.ID,
        'REF': record.REF,
        'ALT': [i.value for i in record.ALT],
        'type': [i.type for i in record.ALT],
        'QUAL': record.QUAL,
        'FILTER': [record.FILTER]
    })
    
    # get VEP annotation values 'CSQ' and create a dataframe
    CSQ_values = record.INFO['CSQ'][0].split("|")
    CSQ_df = pd.DataFrame([CSQ_values], columns=CSQ_fields)
    
    # concatenate both dataframes
    new_row = pd.concat([base_df, CSQ_df[selected_fields]], axis=1)
    
    # append final dataframe with new row
    clinvar_kid_gen2345_vep_anno = pd.concat([clinvar_kid_gen2345_vep_anno, new_row], axis=0)
    
    count = count + 1
    if count % 10000 == 0:
        print(count)
        sys.stdout.flush() 
#     if count > 2000: break

#     if count % 50000 == 0:
#         print(count)
        
# write csv
clinvar_kid_gen2345_vep_anno.to_csv(f"features_clinvar_vars_kid-gen_2345_2024-03-05.csv.gz", index=False, compression='gzip')



