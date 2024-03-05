## TODO: header etc

# basic logic

# export config file
export CONFIG_FILE=/fast/work/users/rankn_c/halbritter/nephro_candidate_score/gene_score/training/config_NCS.yml


# filter Clinvar VCF for kidney-genetic genes with evidence count 2-5
python filter_clinvar_vcf_for_kidney_genes.py > clinvar_vars_kid-gen_2345_2024-03-05.vcf

gzip clinvar_vars_kid-gen_2345_2024-03-05.vcf


# annotate with VEP
# for now with web tool => to be done with command line tool, result: 'clinvar_vars_kid-gen_2345_vep_anno_2024-03-05.vcf.gz'


# extract relevant features from VEP annotation 'CSQ'
python extract_VEP_features.py 
# result: features_clinvar_vars_kid-gen_2345_2024-03-05.csv.gz 




