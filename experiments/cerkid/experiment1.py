#!/usr/bin/env python
# coding: utf-8

# In[3]:
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
import re
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.stats import fisher_exact
from matplotlib.patches import Patch

# import third-party modules
import gseapy as gp

current_date = datetime.now().strftime("%Y-%m-%d")

# TODO: Inheritance pattern is not in long format yet, but only has values: array(['unknown', 'compound_heterozygous_possible_no_pedigree','homozygous']

# ------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------


NGS_predictions_path = (
    "/nephro_candidate_score/gene_score/"
    "predictions/results/"
    "NGS_predictions_ID97_all_2024-03-22.csv.gz"
)

agde_annotation_path = (
    "/data/cephfs-1/work/projects/nephrology-exomes/"
    "multisample-vcf/agde-exomes/results/variantcentrifuge/"
    "agde-exomes_ncs_high_or_lof_or_nmd_rare_mane_select/"
    "agde-exomes.annotated.all.xlsx"
)

kid_gen_path = (
    "/nephro_candidate_score/gene_score/labels/raw/"
    "A_MergeAnalysesSources.2023-10-04.csv.gz"
)

morbid_genes_path = (
    "/nephro_candidate_score/experiments/raw/"
    "MorbidGenes_2025_07.xlsx"
)

# ------------------------------------------------------------------
# NGS predictions
# ------------------------------------------------------------------

# Load NGS prediction table
NGS_all = pd.read_csv(NGS_predictions_path)

# Extract unique genes
all_genes = NGS_all.symbol.unique()


# ------------------------------------------------------------------
# AGDE exomes annotation file
# Filter:
#   ncs_high_or_lof_or_nmd_rare_mane_select
# ------------------------------------------------------------------

# Load annotation table
annot_file1 = pd.read_excel(
    agde_annotation_path,
    engine="openpyxl",
)

# selection 1: ncs_high_or_lof_or_nmd_rare_mane_select
sel1 = pd.read_excel(annot_file1, engine='openpyxl')

# Convert column GT into long format
def convert_GT_into_long_format(df):

    # 1. Split the string by ';' into a list
    df['GT'] = df['GT'].str.split(';')
    
    # 2. Explode the list into individual rows
    df_long = df.explode('GT').reset_index(drop=True)
    
    # 3. Extract the ID and Genotype using Regex
    # Pattern: Capture everything before '(' as patient_id, 
    # and everything inside the '()' as genotype
    df_long[['patient_id', 'genotype']] = df_long['GT'].str.extract(r'^(.*)\((.*)\)$')
    
    # 4. Drop the intermediate helper column
    df_long = df_long.drop(columns=['GT'])
    
    return df_long

# convert GT into long format
sel1_long = convert_GT_into_long_format(sel1)

# filter for N-CS >= 8 
sel1_long_ncs8 = sel1_long.query("nephro_candidate_score >= 8")#.query("control_variant_count <= 2").query("QUAL > 400")


#### Gene counts ####

def get_gene_counts_and_annotate(df):
    """This function groups a df by column 'GENE' and counts variants/Gene. It annotates the genes with kidney-genetics
    evidence count and morbid genes."""

    # get gene counts
    gene_counts = df.groupby('GENE').size().reset_index(name='Count')

    # annotate with kidney-genetis EC
    kid_gen = pd.read_csv(kid_gen_path)
    gene_counts = pd.merge(gene_counts, kid_gen[['approved_symbol', 'evidence_count']], how = 'left', left_on='GENE', right_on='approved_symbol')
    
    # annotate with morbid genes
    morbid = pd.read_excel(morbid_genes_path)
    morbid['morbid_gene'] = 1
    gene_counts = pd.merge(gene_counts, morbid[['Symbol', 'morbid_gene']], how = 'left', left_on='GENE', right_on='Symbol')
    
    # assign color
    conditions = [
        (gene_counts['evidence_count'].between(2, 5)),
        (gene_counts['evidence_count'].between(0, 1)),
        (gene_counts['morbid_gene'] == 1)
    ]
    
    # define corresponding colors
    choices = ['#999999', '#377EB8', '#009E73']
    
    # default color if none of the above
    gene_counts['color'] = np.select(conditions, choices, default='#E69F00' )

    gene_counts = gene_counts.sort_values(by='Count', ascending=False)

    return gene_counts

gene_counts_uncensored = get_gene_counts_and_annotate(sel1_long_ncs8)


# create gene counts plot
df = gene_counts_uncensored.sort_values(by='Count', ascending=False).head(40)

# Plot
plt.figure(figsize=(10, 6))
bars = plt.bar(df['GENE'], df['Count'], color=df['color'])

# Rotate gene labels
plt.xticks(rotation=90)

# Custom legend
legend_elements = [
    Patch(facecolor='#808080', label='Evidence count 2–5'),
    Patch(facecolor='#377EB8', label='Evidence count 0–1'),
    Patch(facecolor='#009E73', label='Morbid gene'),
    Patch(facecolor='#E69F00', label='Other')
]
plt.legend(handles=legend_elements, title='Legend')

# Labels and title
plt.xlabel('Gene')
plt.ylabel('Count')
plt.title('Uncensored Gene Counts Colored by Evidence/Morbid Status')
plt.tight_layout()
plt.show()




