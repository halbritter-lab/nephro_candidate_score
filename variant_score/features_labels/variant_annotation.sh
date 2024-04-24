## TODO: header etc

# basic logic

# execution from nephro_candidate_score/variant_score

# export config file
export CONFIG_FILE=/fast/work/users/rankn_c/halbritter/nephro_candidate_score/gene_score/training/config_NCS.yml


## download ClinVar VCF file for GRCh38
# URL of the ClinVar VCF file for GRCh38
clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20240301.vcf.gz"

# Destination file path
destination_path = "raw/clinvar_20240301.vcf.gz"

# Download the file
# urllib.request.urlretrieve(clinvar_url, destination_path)

# Check if the download was successful
if os.path.exists(destination_path):
    print("ClinVar VCF file downloaded successfully.")
else:
    print("Failed to download ClinVar VCF file.")



# filter Clinvar VCF for kidney-genetic genes with evidence count 2-5
python features_labels/filter_clinvar_vcf_for_kidney_genes.py > features_labels/clinvar_vars_kid-gen_2345_2024-03-07.vcf

gzip features_labels/clinvar_vars_kid-gen_2345_2024-03-07.vcf


# annotate with VEP
# for now with web tool => to be done with command line tool, result: 'clinvar_vars_kid-gen_2345_vep_anno_2024-03-05.vcf.gz'


# extract relevant features from VEP annotated VCF
python extract_variant_features.py 
# results: 
# - raw_features_clinvar_vars_kid-gen_2345_prefiltered_{datetime.today().strftime('%Y-%m-%d')}.csv.gz" # only filtered for kid-gen2345 genes in 'filter_clinvar_vcf_for_kidney_genes.py' => but through VEP annotation still vars in genes other than kid-gen2345 genes
# - raw_features_clinvar_vars_kid-gen_2345_{datetime.today().strftime('%Y-%m-%d')}.csv.gz # again filtered for kid-gen2345 genes => only variants in kid-gen2345 genes



