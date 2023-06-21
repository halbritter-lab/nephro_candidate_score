# Nephro Candidate Score - a tool for automatic prioritization of variants in kidney disease patients

Welcome to the GitHub repository  "Nephro Candidate Score". The Nephro Candidate Score (NCS) is a tool for standardized and automatic prioritization of candidate variants in kidney disease patients. It integrates information from multiple open-source databases and scores variants based on a machine learning algorithm.

## Table of contents

- [Overview and Methods](#overview-and-methods)
- [Usage](#usage)
- [Documentation](#documentation)
- [File structure](#file-structure)
- [License](#license)
- [Creators](#creators-and-contributors)
- [Contact](#contact)


## Overview and Methods
TODO


## Usage
A webtool for single variant scoring is available on TODO.
A command line tool for scoring .vcf files is available on TODO.

## Documentation
TODO


## File Structure

The repository has the following structure:

```
.
├── gene_score/
│   ├── 01_PanelApp/
│   │   ├── data/
│   │   ├── results/
│   │   └── 01_PanelApp.R
│   ├── 02_Literature/
│   │   ├── data/
│   │   ├── results/
│   │   └── 02_Literature.R
│   ├── 03_DiagnosticPanels/
│   │   ├── data/
│   │   ├── results/
│   │   └── 03_DiagnosticPanels.R
│   ├── 04_HPO/
│   │   ├── data/
│   │   ├── results/
│   │   └── 04_HPO.R
│   ├── 05_PubTator/
│   │   ├── data/
│   │   ├── results/
│   │   └── 05_PubTator.R
│   ├── MergeAnalysesSources.R
│   └── AnnotateMergedTable.R
└── functions/
    ├── blueprintgenetics-functions.R
    ├── hgnc-functions.R
    ├── hpo-functions.R
    ├── natera-functions.R
    ├── NCBI-datasets-v2-API-functions.R
    ├── phantomjs-functions.R
    └── PubTator-functions.R
```

- The `analyses/` directory contains the R scripts for different analyses.
- The `functions/` directory contains the necessary functions for HGNC processing.
- The `data/` sub-directory in each analysis folder stores the input data files, including the publication-specific files and the curated overview Excel table.
- The `results/` sub-directory in each analysis folder stores the generated results.


## License

This project is licensed under the terms of the MIT license. For more information, please refer to the [License](LICENSE.md) file.


## Creators and contributors
**Nina Rank**

- <https://github.com/ninarank>
- <https://orcid.org/0000-0002-5984-4836>


**Bernt Popp**

- <https://twitter.com/berntpopp>
- <https://github.com/berntpopp>
- <https://orcid.org/0000-0002-3679-1081>
- <https://scholar.google.com/citations?user=Uvhu3t0AAAAJ>

**Constantin Aaron Wolff**

- <https://github.com/ConstantinWolff>
- <https://orcid.org/0000-0002-3277-4559>

**Jan Halbritter**

- <https://orcid.org/0000-0002-1377-9880>
- <https://scholar.google.com/citations?user=Jt1S5fkAAAAJ>

## Contact

If you have any questions, suggestions, or feedback, please feel free to [contact us](contact.md).
