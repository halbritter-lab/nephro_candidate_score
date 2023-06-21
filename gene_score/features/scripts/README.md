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

#### Number of extracted features
TODO
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

#### Number of extracted features
2
---

### gnomad.R
#### Datasource: 
gnomAD (https://gnomad.broadinstitute.org/, https://doi.org/10.1038/s41586-020-2308-7)

#### Description
This script extracts gnomAD gene constraint metrics, which are published here: https://gnomad.broadinstitute.org/downloads.
  
#### Required libraries
- `tidyverse`: for data processing.
- `utils`: for downloading data.

#### Number of extracted features
TODO


