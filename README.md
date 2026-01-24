# `NucFlag`
[![CI](https://github.com/logsdon-lab/NucFlag/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/NucFlag/actions/workflows/main.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/nucflag)](https://pypi.org/project/nucflag/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-blue.svg?style=flat)](http://bioconda.github.io/recipes/nucflag/README.html)

Generates nucleotide frequency plots and genome misassembly BED files. Fork of [`NucFreq`](https://github.com/mrvollger/NucFreq).

![Labeled Misassemblies](docs/imgs/misassemblies.png)

## Quickstart
```bash
pip install nucflag
# Or conda install nucflag
```

> [!NOTE]
> NucFlag v1.0.0-alpha.1 is not installable via `bioconda` at the moment.

Align long-reads to assembly.
```bash
# Align long reads to assembly generated from those reads.
minimap2 -x lr:hqae -I 8G asm.fa.gz reads.fq.gz | samtools view -bh -o asm.bam
samtools index asm.bam
```

Detect misassemblies.
```bash
# Call putative misassemblies from read alignments on a whole genome.
nucflag call -i asm.bam -f asm.fa.gz -o misassemblies.bed
# Or on a set of regions...
nucflag call -i asm.bam -f asm.fa.gz -b regions.bed -o misassemblies.bed
# Also runs on ONT read alignments.
nucflag call -i asm_ont_r9.bam -f asm.fa.gz -x ont_r9 -o misassemblies.bed
nucflag call -i asm_ont_r10.bam -f asm.fa.gz -x ont_r10 -o misassemblies.bed
# Provide a configfile for finer control.
nucflag call -i asm.bam -f asm.fa.gz -o misassemblies.bed -c config.toml
```

Visualize misassemblies in a number of ways.
```bash
# Generate NucFreq plots.
nucflag call -i asm.bam -f asm.fa.gz -d plots
# And add any number of tracks...
nucflag call -i asm.bam -f asm.fa.gz -d plots --tracks repeatmasker.bed segdups.bed
# Or generate bigWigs of specific signals and then merge them with `bigtools`.
# For use in IGV or other genome browsers.
nucflag call -i asm.bam -f asm.fa.gz -d plots \
    --output_pileup_dir bigwigs \
    --add_pileup_data cov mismatch mapq
bigwigmerge -l <(find bigwigs -name "*_first.bw") merged_first.bw
```

Generate status BED or breakdown showing distribution of assembly issues.
```bash
nucflag status -i misassemblies.bed > status.bed
# Or plot by "length"
nucflag breakdown -i misassemblies.bed -o breakdown -t percent
```

Estimate QV from BED file.
```bash
nucflag qv -i misassemblies.bed > qv.bed
```

Generate ideogram.
```bash
nucflag ideogram -i misassemblies.bed -o ideogram
# Add cytobands with a BED file with chrom, start, end, name, and btype.
nucflag ideogram -i misassemblies.bed -c cytobands.bed -o ideogram
```

Get consensus misassembly calls by intersection.
```bash
nucflag consensus -i nucflag_ont.bed nucflag_hifi.bed hmm_flagger_hifi.bed hmm_flagger_ont.bed > consensus.bed
```

Generate sample config.
```bash
nucflag config -x ont_r10 > config_r10.toml
```

## Input
* BAM file of PacBio HiFi, ONT R9, or ONT R10 reads aligned to an assembly.
* (Recommended) Assembly.
* (Optional) BED file of regions.

## Output
* BED file of misassemblies.
* (Optional) Plots with coverage and mismatch pileup with misassemblies flagged.
* (Optional) bigWigs of pileup signals.
* (Optional) BED file of assembly status.

## [Documentation](https://github.com/logsdon-lab/NucFlag/wiki)
Read the docs at the `NucFlag` [wiki](https://github.com/logsdon-lab/NucFlag/wiki) for more information.
