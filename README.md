# `NucFlag`
[![CI](https://github.com/logsdon-lab/NucFlag/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/NucFlag/actions/workflows/main.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/nucflag)](https://pypi.org/project/nucflag/)

Fork of [`NucFreq`](https://github.com/mrvollger/NucFreq). Script for making nucleotide frequency plots and marking misassemblies.

![Labeled Misassemblies](docs/imgs/misassemblies.png)

## Usage
```bash
pip install nucflag
```

```
usage: nucflag [-h] -i INFILE [-b INPUT_REGIONS] [-d OUTPUT_PLOT_DIR] [--output_cov_dir OUTPUT_COV_DIR] [-o OUTPUT_MISASM] [-s OUTPUT_STATUS] [-t THREADS] [-p PROCESSES] [-c CONFIG] [--ignore_regions IGNORE_REGIONS] [--overlay_regions [OVERLAY_REGIONS ...]]

Use per-base read coverage to classify/plot misassemblies.

options:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input bam file or per-base coverage tsv file with 3-columns (position, first, second). If a bam file is provided, it must be indexed. (default:
                        None)
  -b INPUT_REGIONS, --input_regions INPUT_REGIONS
                        Bed file with regions to check. (default: None)
  -d OUTPUT_PLOT_DIR, --output_plot_dir OUTPUT_PLOT_DIR
                        Output plot dir. (default: None)
  --output_cov_dir OUTPUT_COV_DIR
                        Output coverage dir. Generates gzipped coverage bed files per region. (default: None)
  -o OUTPUT_MISASM, --output_misasm OUTPUT_MISASM
                        Output bed file with misassembled regions. (default: <_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>)
  -s OUTPUT_STATUS, --output_status OUTPUT_STATUS
                        Bed file with status of contigs. With format: contig start end misassembled|good (default: None)
  -t THREADS, --threads THREADS
                        Threads for reading bam file. (default: 4)
  -p PROCESSES, --processes PROCESSES
                        Processes for classifying/plotting. (default: 4)
  -c CONFIG, --config CONFIG
                        Additional threshold/params as toml file. (default: {'first': {'thr_min_peak_horizontal_distance': 1, 'thr_min_peak_width': 20,
                        'thr_min_valley_horizontal_distance': 1, 'thr_min_valley_width': 3, 'thr_peak_height_std_above': 4, 'thr_valley_height_std_below': 2,
                        'thr_misjoin_valley': 3, 'valley_group_distance': 500, 'peak_group_distance': 500}, 'second': {'thr_min_perc_first': 0.1,
                        'thr_peak_height_std_above': 3, 'group_distance': 30000, 'thr_min_group_size': 3, 'thr_collapse_het_ratio': 0.2}, 'gaps':
                        {'thr_max_allowed_gap_size': 0}})
  --ignore_regions IGNORE_REGIONS
                        Bed file with regions to ignore. With format: contig|all start end absolute|relative (default: None)
  --overlay_regions [OVERLAY_REGIONS ...]
                        Overlay additional regions as 4-column bedfile alongside coverage plot. (default: None)
```

### Input
A BAM file of an alignment of PacBio HiFi reads to an assembly.

> [!IMPORTANT]
> All assembly contigs, including contaminants, should be included. Omission of these contigs will cause misalignments of reads to elsewhere in the assembly.

Secondary and partial alignments should be removed using SAMtools flag 2308.

### Configuration

#### `--config`
Configuration can be provided in the form of a `toml` file and the `--config` flag.

```bash
nucflag -i test/HG00096_hifi_test.bam -b test/test.bed -c config.toml
```

```toml
[first]
# Min horizontal distance between peaks.
thr_min_peak_horizontal_distance = 1
# Min width of peak to consider.
thr_min_peak_width = 20
# Min horizontal distance between valleys.
thr_min_valley_horizontal_distance = 1
# Min width of valley to consider.
thr_min_valley_width = 10
# Number of std above mean to include peak.
thr_peak_height_std_above = 4
# Number of std below mean to include valley.
thr_valley_height_std_below = 3
# Valleys with coverage below this threshold are considered misjoins.
# If float:
# * ex. 0.1 => Valley where min is less than or equal to 10% of mean
# If int:
# * ex. 2 => Valleys with min below 2
thr_misjoin_valley = 0.1
# Group consecutive positions allowing a maximum gap of x.
# Larger value groups more positions.
valley_group_distance = 500
peak_group_distance = 500

[second]
# Percent threshold of most freq base to allow second most freq base
# 10 * 0.1 = 1 so above 1 is allowed.
thr_min_perc_first = 0.1
# Number of std above mean to include peak.
thr_peak_height_std_above = 3
# Group consecutive positions allowing a maximum gap of x.
# Larger value groups more positions.
group_distance = 30_000
# Min group size.
thr_min_group_size = 5
# Het ratio to consider second group a collapse if no overlaps in peaks found.
thr_collapse_het_ratio = 0.1

[gaps]
# Allow gaps up to this length.
thr_max_allowed_gap_size = 1000
```

#### `--ignore_regions`

Specific regions can be ignored via a 5-column BED file and the `--ignore_regions` flag.

```bash
nucflag -i NA20847.bam -b region.bed --ignore_regions ignore.bed
```

|contig|start|stop|desc|action|
|-|-|-|-|-|
|all|0|500000|-|ignore:relative|
|all|-500000|0|-|ignore:relative|
|ctg_name|0|500000|-|ignore:absolute|

Two options are posible for `ignore` action.
* `relative`
  * Coordinates will be relative to the full length of the contig and a position of `0` will anchor which side to ignore.

* `absolute`
  * Misassemblies within the bounds of the contig and the coordinates will be ignored.

> Ex. Ignores misasemblies within 500kbp of the edges of all contigs.
```
all	0	500000	pericentromere	ignore:relative
all	-500000	0	pericentromere	ignore:relative
```

> Ex. Ignores misasemblies between 10bp to 50bp on `hap1-0000001`.
```
hap1-0000001	10	50	test	ignore:absolute
```

#### `--overlap_regions`
Regions can also be added as tracks via a 5-column BED file and the `--overlap_regions` flag.
```bash
nucflag -i NA20847.bam -b region.bed --overlap_regions NA20847_repeatmasker.bed
```
```
NA20847_rc-chr3_haplotype2-0000105	89882645	89883610	LTR	plot
NA20847_rc-chr3_haplotype2-0000105	89883812	89884058	LTR	plot
NA20847_rc-chr3_haplotype2-0000105	89884065	89884093	Simple_repeat	plot
NA20847_rc-chr3_haplotype2-0000105	89884094	89884130	Simple_repeat	plot
...
```

![Overlap Bed](docs/imgs/overlap.png)


## Workflow
For an end-to-end workflow, see [`Snakemake-NucFlag`](https://github.com/logsdon-lab/Snakemake-NucFlag).

## Build
To build from source.
```bash
git clone git@github.com:logsdon-lab/NucFlag.git && cd NucFlag
make venv && make build && make install
```

## Test
Test BAM filtered from merged alignment of:
* PacBio HiFi reads from HGSVC sample `HG00096`.
  * https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20220831_JAX_HiFi/HG00096/
* Verkko v1.4.1 combined assembly for HGSVC sample `HG00096`
  * https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20240201_verkko_batch3/assemblies/HG00096/

To run tests:
```bash
source venv/bin/activate
make test
```

Or try the test example directly.
```bash
nucflag -i test/standard/HG00096_hifi.bam -b test/standard/region.bed -c test/config.toml
```
```
haplotype2-0000133:3021508-8691473      3314093 3324276 MISJOIN
haplotype2-0000133:3021508-8691473      4126277 4142368 COLLAPSE_VAR
haplotype2-0000133:3021508-8691473      4566798 4683011 GAP
haplotype2-0000133:3021508-8691473      5737835 5747246 MISJOIN
haplotype2-0000133:3021508-8691473      6067838 6072601 COLLAPSE_VAR
haplotype2-0000133:3021508-8691473      6607947 6639102 MISJOIN
haplotype2-0000133:3021508-8691473      7997560 8069465 COLLAPSE_VAR
```

Test workflow using `data/` dir. Requires:
1. `bam` files of alignment of HGSVC3 HiFi reads to full assembly.
2. `bed` files of centromere + 500kbp regions.

```bash
snakemake \
-s test/workflow/Snakefile \
-j 12 \
--executor cluster-generic \
--cluster-generic-submit-cmd "bsub -q epistasis_normal -n {threads} -o /dev/null" \
--use-conda -p
```

## Cite
- **Vollger MR**, Dishuck PC, Sorensen M, Welch AE, Dang V, Dougherty ML, et al. Long-read sequence and assembly of segmental duplications. Nat Methods. 2019;16: 88–94. doi:10.1038/s41592-018-0236-3
- **Mc Cartney AM**, Shafin K, Alonge M, Bzikadze AV, Formenti G, Fungtammasan A, et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. bioRxiv. 2021. p. 2021.07.02.450803. doi:10.1101/2021.07.02.450803
  * Citing `hetDetection.R`
