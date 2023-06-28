# NCS - Label aquisition

The scripts in this directory retrieve the labels (positive and dispensable genes) that are used for training the machine learning algorithm in a later step.


## Script description
### dispensable_genes.R
#### Datasources: 
- gnomAD     
(Karczewski KJ et al. Genome Aggregation Database Consortium; Neale BM, Daly MJ, MacArthur DG. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature (2020). doi: 10.1038/s41586-020-2308-7), Supplementary Dataset 7  
- OMIM (https://www.omim.org/)   
- kidney-genetics (https://github.com/halbritter-lab/kidney-genetics/)

#### Description
This scripts determines dispensable genes - genes with homozygous loss of function variants that are not associated with any disease. It first downloads a list of homozygous knockout genes from the gnomAD publication (s. above). Additionally, it downloads the genemap2.txt file from OMIM - a tab-delimited file containing OMIM's Synopsis of the Human Gene Map including additional information such as genomic coordinates and inheritance. From the list of homozygous ko genes, it removes those which have a 'Phenotype' entry in OMIM or an entry in 'kidney-genetics'. 'Kidney-genetics' is a reproducible and curated database of kidney disease related genes. 

#### Requirements
- `tidyverse`: for data processing.
- `readr`: for sending GET and POST requests.
- `jsonlite`: for parsing .json files.   

This script sources an R-script from the github repository kidney-genetics.   
Downloading the genemap2.txt file from OMIM requires registration at the OMIM website (https://www.omim.org/downloads). The personal download link has to be added to the config file of this repository.


---

### positive_genes.R
#### Datasources: 
- kidney-genetics (https://github.com/halbritter-lab/kidney-genetics/)

#### Description
This scripts sources genes from the github repository 'kidney-genetics'. 'Kidney-genetics' is a reproducible and curated database of kidney disease related genes. 
 
#### Requirements
- `tidyverse`: for data processing.
- `utils`: for downloading files.   
---








