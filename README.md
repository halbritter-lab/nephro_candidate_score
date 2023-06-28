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
│   ├── features/
│   │   ├── raw/
│   │   ├── results/
│   │   ├── scripts/
│   │   │   ├── cellxgene.R
│   │   │   ├── exon_and_prom_conservation.R
│   │   │   ├── gnomad.R
│   │   │   ├── gtex.R
│   │   │   ├── kidney_network.R
│   │   │   ├── mgi_mpo.R
│   │   │   ├── nephrogenesis_atlas.R
│   │   │   ├── paralogues.R
│   │   │   ├── promoter_CpG_o2e_ratio.R
│   ├── labels/
│   │   ├── raw/
│   │   ├── results/
│   │   ├── scripts/
│   │   │   ├── dispensable_genes.R
│   │   │   ├── positive_genes.R
└── variant_score/
```

- The `gene_score/` directory contains TODO
- The `functions/` directory contains TODO
- The `data/` TODO
- The `results/` TODO


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
