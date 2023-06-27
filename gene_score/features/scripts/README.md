# NCS - Feature aquisition

The scripts in this directory retrieve features from different open-source databases that are used for training the machine learning algorithm in a later step.


## Script description
### cellxgene.R
#### Datasource: 
CZ CELLxGENE Discover (https://cellxgene.cziscience.com/, Chan Zuckerberg Initiative. (n.d.).)

#### Description
This script extracts ssRNA expression values from cellxgene via the cellxgene API. In specific, it pulls the expression data from TODO

#### Required libraries
- `tidyverse`: for data processing.
- `httr`: for sending GET and POST requests.
- `jsonlite`: for parsing .json files.

#### Number of extracted features: TODO
---

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

#### Number of extracted features: 2
---

### gnomad.R
#### Datasource: 
gnomAD (https://gnomad.broadinstitute.org/, https://doi.org/10.1038/s41586-020-2308-7)

#### Description
This script extracts gnomAD gene constraint metrics, which are published here: https://gnomad.broadinstitute.org/downloads.
  
#### Required libraries
- `tidyverse`: for data processing.
- `utils`: for downloading data.

#### Number of extracted features: TODO

---

### gtex.R
#### Datasource: 
The Human Protein Atlas (https://www.proteinatlas.org/downloads, RNA GTEx tissue gene data). The data was originally obtained from GTEx and is based on The Human Protein Atlas version 23.0 and Ensembl version 109.

#### Description
This script extracts bulk RNA tissue normalized expression ("nTPM") values from 'The Human Protein Atlas' (primary data source GTEx). It aggregates the nTPM values for the different brain regions to one brain_nTPM_med value. Additionally, it calculates the normalized &tau;  tissue specificity index for each gene according to Yanai et al (Yanai et al. Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification, Bioinformatics, Volume 21, Issue 5, March 2005, Pages 650?659, https://doi.org/10.1093/bioinformatics/bti042). &tau; is a metric to quantify the tissue specificity of gene expression. It ranges from 0 to 1, with values closer to 1 indicating higher tissue specificity. &tau; = 1 suggests that the gene is exclusively expressed in a single tissue, while a &tau; closer to 0 indicates a more widespread or ubiquitous expression across multiple tissues.

#### Required libraries
- `tidyverse`: for data processing.
- `utils`: for downloading data.

#### Number of extracted features: 29 + 1
---

### kidney_network.R
#### Datasource: 
KidneyNetwork (https://kidney.genenetwork.nl)

#### Description
This script extracts gene-pathway prediction z-scores of bonferroni corrected significantly predicted pathways from KidneyNetwork. It sums up the positive (negative z-scores) of all children of each of the following HPO-terms:
1. "Abnormal renal morphology", HP:0012210  
2. "Abnormal renal physiology", HP:0012211  
3. "Abnormality of the urinary system", HP:0000079 (includes the upper two HPO-terms)

#### Required libraries
- `tidyverse`: for data processing.
- `utils`: for downloading data.

#### Number of extracted features: 6


---

### mgi_mpo.R
#### Datasource: 
Mouse Genome Informatics (https://www.informatics.jax.org)

#### Description
This script identifies which mouse genotypes (homozygous or heterozygous knock-out mice) are  associated with mouse phenotype ontology (MPO) term "MP:0005367" (= "renal/urinary system phenotype") or any of its children. It then annotates the respective mouse genes with the corresponding human ortholog.

#### Required libraries
- `tidyverse`: for data processing.
- `utils`: for downloading data.
- `jsonlite`: for parsing .json files.

#### Number of extracted features: 2 (heterozygous/homozygous)





