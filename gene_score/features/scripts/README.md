# NCS - Feature aquisition

The scripts in this directory retrieve features from different open-source databases that are used for training the machine learning algorithm in a later step.


## Script description
### cellxgene.R
#### Datasource: 
CZ CELLxGENE Discover (https://cellxgene.cziscience.com/, Chan Zuckerberg Initiative. (n.d.).)

#### Description
This script extracts ssRNA expression values from cellxgene via the cellxgene API.

#### Required libraries
- `tidyverse`: for data processing.
- `httr`: for sending GET and POST requests.
- `jsonlite`: for parsing .json files.

#### Number of extracted features
TODO





