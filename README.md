# Nephro Candidate Score - a tool for automatic prioritization of variants in kidney disease patients

Welcome to the GitHub repository  "Nephro Candidate Score". The Nephro Candidate Score (NCS) is a tool for standardized and automatic prioritization of candidate variants in kidney disease patients. It integrates information from multiple open-source databases and scores variants based on a machine learning algorithm.

## Table of contents

- [Overview and Methods](#overview-and-methods)
- [Usage](#usage)
- [File structure](#file-structure)
- [License](#license)
- [Creators](#creators-and-contributors)
- [Contact](#contact)


## Overview and Methods
TODO


## Usage
The NCS is publicly available and accessible. A web tool version for single variant scoring is available on TBA.
A command line tool for scoring .vcf files is available on TBA.


## File Structure

The repository has the following structure:


```
.
├── gene_score/
│   ├── features/
│   │   ├── raw/
│   │   ├── results/
│   │   ├── scripts/
│   │   │   ├── cellxgene.R
│   │   │   ├── descartes.R
│   │   │   ├── exon_CpG_o2e_ratio.R
│   │   │   ├── exon_and_prom_conservation.R
│   │   │   ├── gnomad.R
│   │   │   ├── gtex.R
│   │   │   ├── kidney_network.R
│   │   │   ├── mgi_mpo.R
│   │   │   ├── nephrogenesis_atlas.R
│   │   │   ├── paralogues.R
│   │   │   ├── promoter_CpG_o2e_ratio.R
│   │   ├── helper_functions.R
│   ├── labels/
│   │   ├── raw/
│   │   ├── results/
│   │   ├── scripts/
│   │   │   ├── dispensable_genes.R
│   │   │   ├── positive_genes.R
│   ├── raw/
│   ├── training/
│   ├── gene_score_preprocessing_master.R
│   ├── helper_functions.R
│   ├── hgnc_functions.R
└── variant_score/
```

- The `gene_score/` directory contains all scripts relevant for data acquisition and machine learning development of the gene score. Due to large file sizes, the raw data is not completely stored on GitHub.
- The `variant_score/` TBA


## License

This project is licensed under the terms of the MIT license. For more information, please refer to the [License](LICENSE.md) TBA file.


## Creators and contributors
**Nina Rank**

- <https://github.com/ninarank>
- <https://orcid.org/0000-0002-5984-4836>


**Bernt Popp**

- <https://twitter.com/berntpopp>
- <https://github.com/berntpopp>
- <https://orcid.org/0000-0002-3679-1081>
- <https://scholar.google.com/citations?user=Uvhu3t0AAAAJ>


**Soeren Lukassen**
- <https://github.com/slukassen>
- <https://orcid.org/0000-0001-7045-6327>
- <https://scholar.google.com/citations?user=wWiHGZkAAAAJ&hl=en>


**Constantin Aaron Wolff**

- <https://github.com/ConstantinWolff>
- <https://orcid.org/0000-0002-3277-4559>

**Jan Halbritter**

- <https://orcid.org/0000-0002-1377-9880>
- <https://scholar.google.com/citations?user=Jt1S5fkAAAAJ>

## Contact

If you have any questions, suggestions, or feedback, please feel free to contact nina.rank@charite.de.
