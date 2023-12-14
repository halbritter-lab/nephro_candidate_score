# NCS - Feature aquisition

The scripts in this directory retrieve features from different open-source databases that are used for training the machine learning algorithm in a later step.


## Script description
### cellxgene.R
#### Datasource: 
CZ CELLxGENE Discover (https://cellxgene.cziscience.com/, Chan Zuckerberg Initiative. (n.d.).)

#### Description
This script extracts ssRNA expression values from cellxgene via the cellxgene API. In specific, it pulls the expression data from the following dataset: 
https://cellxgene.dev.single-cell.czi.technology/collections/a98b828a-622a-483a-80e0-15703678befd by Markus Bitzer, University of Michigan.
The scripts extracts 'me' (Gene expression) and 'pc' (% Expressed in Cells) values of all kidney cell types. 

#### Required libraries
- `tidyverse`: for data processing.
- `httr`: for sending GET and POST requests.
- `jsonlite`: for parsing .json files.
- `progress`: for progress bars.
- `R.utils`: for gzipping files.
- `config`: for loading configurations.


#### Number of extracted features: 26
---

### exon_CpG_o2e_ratio.R
#### Datasource: 
Ensembl Biomart (doi:10.1093/nar/gkab1049)

#### Description
This script extracts positional information from canonical transcripts of protein coding genes from Ensembl Biomart. It calculates the observed-to-expected-CpG-ratio of each exon of each gene. It the calculates the average exon observed-to-expected-CpG-ratio per gene. 

#### Required libraries
- `tidyverse`: for data processing.
- `biomaRt`: for extracting data from Ensembl biomart.
- `progress`: for progress bars.
- `BSgenome.Hsapiens.UCSC.hg38`: for extracting the exon sequences. 
- `R.utils`: for gzipping files.
- `config`: for loading configurations.

#### Number of extracted features: 1


### exon_and_prom_conservation.R
#### Datasource: 
Ensembl Biomart (doi:10.1093/nar/gkab1049)

#### Description
This script extracts exon information from Ensembl Biomart and calculates the average conservation PhastCons scores of all coding exons and of the promoter region of each canonical transcript, respectively. The canonical transcript is defined as the Ensembl canonical transcript (https://www.ensembl.org/info/genome/genebuild/canonical.html). 
The promoter region is defined as +-2kb around the transcription start site. The PhastCons score ranges from 0 to 1 and represents the probability that a given nucleotide is conserved.  
  
#### Required libraries
- `tidyverse`: for data processing.
- `biomaRt`: for extracting data from Ensembl biomart.
- `progress`: for progress bars for long calculation time.
- `phastCons100way.UCSC.hg38`: contains the PhastCons scores. 
- `R.utils`: for gzipping files.
- `config`: for loading configurations.


#### Number of extracted features: 2
---

### gnomad.R
#### Datasource: 
gnomAD, v2.1.1 (https://gnomad.broadinstitute.org/, https://doi.org/10.1038/s41586-020-2308-7)

#### Description
This script extracts gnomAD gene constraint metrics, which are published here: https://gnomad.broadinstitute.org/downloads.
  
#### Required libraries
- `tidyverse`: for data processing.
- `utils`: for downloading data.
- `R.utils`: for gzipping files.
- `config`: for loading configurations.

#### Number of extracted features: 34
---

### gtex.R
#### Datasource: 
The Human Protein Atlas (https://www.proteinatlas.org/downloads, RNA GTEx tissue gene data). The data was originally obtained from GTEx and is based on The Human Protein Atlas version 23.0 and Ensembl version 109.

#### Description
This script extracts bulk RNA tissue normalized expression ("nTPM") values from 'The Human Protein Atlas' (primary data source GTEx). It aggregates the nTPM values for the different brain regions to one brain_nTPM_med value. Additionally, it calculates the normalized &tau;  tissue specificity index for each gene according to Yanai et al (Yanai et al. Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification, Bioinformatics, Volume 21, Issue 5, March 2005, Pages 650?659, https://doi.org/10.1093/bioinformatics/bti042). &tau; is a metric to quantify the tissue specificity of gene expression. It ranges from 0 to 1, with values closer to 1 indicating higher tissue specificity. &tau; = 1 suggests that the gene is exclusively expressed in a single tissue, while a &tau; closer to 0 indicates a more widespread or ubiquitous expression across multiple tissues.

#### Required libraries
- `tidyverse`: for data processing.
- `utils`: for downloading data.
- `R.utils`: for gzipping files.
- `config`: for loading configurations.


#### Number of extracted features: 29 + 1
---

### mgi_mpo.R
#### Datasource: 
Mouse Genome Informatics (https://www.informatics.jax.org)

#### Description
This script identifies which mouse genotypes (homozygous or heterozygous knock-out mice) are  associated with mouse phenotype ontology (MPO) term "MP:0005367" (= "renal/urinary system phenotype") or any of its children. It then annotates the respective mouse genes with the corresponding human ortholog. The resulting feature has 4 values:
2: heterozygous knock-out mice of this gene are associated with mouse phenotype ontology (MPO) term "MP:0005367" or any of its children
1: homozygous knock-out mice of this gene are associated with mouse phenotype ontology (MPO) term "MP:0005367" or any of its children
0: no phenotypes at all exist for this gene
-1: phenotypes exist for this gene, but not mouse phenotype ontology (MPO) term "MP:0005367" or any of its children 


#### Required libraries
- `tidyverse`: for data processing.
- `utils`: for downloading data.
- `jsonlite`: for parsing .json files.
- `R.utils`: for gzipping files.
- `config`: for loading configurations.

#### Number of extracted features: 1

---

### paralogues.R
#### Datasource: 
Ensembl Biomart (doi:10.1093/nar/gkab1049)

#### Description
This script extracts paralogues of protein-coding genes from Ensembl Biomart. It determines the number of close paralogues above the 95th, the 85th and the 75th percentile (Target \%ID and Query \%ID).   
Target \%ID = percentage of identical amino acids in the paralogue compared with the gene of interest    
Query \%ID = percentage of identical amino acids in the gene of interest compared with the paralogue   

#### Required libraries
- `tidyverse`: for data processing.
- `biomaRt`: for extracting data from Ensembl biomart.
- `R.utils`: for gzipping files.
- `config`: for loading configurations.

#### Number of extracted features: 3

---
### promoter_CpG_o2e_ratio.R
#### Datasource: 
Ensembl Biomart (doi:10.1093/nar/gkab1049)

#### Description
This script extracts positional information from canonical transcripts of protein coding genes from Ensembl Biomart. It the calculates the observed-to-expected-CpG-ratio of the promoter region of each gene. The promoter region is defined as +-2kb around the transcription start site. 

#### Required libraries
- `tidyverse`: for data processing.
- `biomaRt`: for extracting data from Ensembl biomart.
- `BSgenome.Hsapiens.UCSC.hg38`: for extracting the promoter sequences. 
- `R.utils`: for gzipping files.
- `config`: for loading configurations.

#### Number of extracted features: 1






