# `NucFlag`
[![CI](https://github.com/logsdon-lab/NucFlag/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/NucFlag/actions/workflows/main.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/nucflag)](https://pypi.org/project/nucflag/)

Generates nucleotide frequency plots and genome misassembly BED files. Fork of [`NucFreq`](https://github.com/mrvollger/NucFreq).

![Labeled Misassemblies](docs/imgs/misassemblies.png)

## Input
* BAM file of PacBio HiFi reads to an assembly.
* (Optional) BED file of regions.

## Output
* BED file of misassemblies.
* Plot with coverage of 1st and 2nd most common bases and misassemblies flagged.

## [Documentation](https://github.com/logsdon-lab/NucFlag/wiki)
Read the docs at the `NucFlag` [wiki](https://github.com/logsdon-lab/NucFlag/wiki) for more information.
