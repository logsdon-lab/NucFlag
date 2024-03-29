# NucFreq Misassembly Classifier
[![CI](https://github.com/logsdon-lab/NucFreq/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/NucFreq/actions/workflows/main.yml)

Fork of [`NucFreq`](https://github.com/mrvollger/NucFreq). Script for making nucleotide frequency plots and marking misassemblies.

![Labeled Misassemblies](docs/imgs/misassemblies.png)

## Usage
Install from GitHub.
```bash
pip install git+https://github.com/logsdon-lab/NucFreq.git
```

```
usage: NucPlot.py [-h] [-i INPUT_BAM] [-b INPUT_BED] [-d OUTPUT_PLOT_DIR] [-o OUTPUT_BED] [-r [REGIONS ...]] [-t THREADS] [-p PROCESSES]

Use per-base read coverage to classify/plot misassemblies.

options:
  -h, --help            show this help message and exit
  -i INPUT_BAM, --input_bam INPUT_BAM
                        Input bam file. Must be indexed. (default: None)
  -b INPUT_BED, --input_bed INPUT_BED
                        Bed file with regions to plot. (default: None)
  -d OUTPUT_PLOT_DIR, --output_plot_dir OUTPUT_PLOT_DIR
                        Output plot dir. (default: None)
  -o OUTPUT_BED, --output_bed OUTPUT_BED
                        Output bed file with misassembled regions. (default: <_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>)
  -r [REGIONS ...], --regions [REGIONS ...]
                        Regions with the format: (.*):(\d+)-(\d+) (default: None)
  -t THREADS, --threads THREADS
                        Threads for reading bam file. (default: 4)
  -p PROCESSES, --processes PROCESSES
                        Processes for classifying/plotting. (default: 4)
```

### Configuration
Configuration can be provided in the form of a `toml` file.

```bash
nucfreq -i test/HG00096_hifi_test.bam -b test/test.bed -c config.toml
```

```toml
[first]
# Bases to add to region bounds
added_region_bounds = 0
# Min horizontal distance between peaks.
thr_min_peak_horizontal_distance = 100_000
# Min width of peak to consider.
thr_min_peak_width = 20
# Min horizontal distance between valleys.
thr_min_valley_horizontal_distance = 100_000
# Min width of valley to consider.
thr_min_valley_width = 10
# Number of std above mean to include peak.
thr_peak_height_std_above = 3.5
# Number of std below mean to include valley.
thr_valley_height_std_below = 3

[second]
# Percent threshold of most freq base to allow second most freq base
# 10 * 0.1 = 1 so above 1 is allowed.
thr_min_perc_first = 0.07
# Number of std above mean to include peak.
thr_peak_height_std_above = 3
# Group consecutive positions allowing a maximum gap of x.
# Larger value groups more positions.
group_distance = 25_000
# Min group size.
thr_min_group_size = 5
# Min group len from starting position to ending position.
thr_min_group_len = 20_000
# Het ratio to consider second group a collapse if no overlaps in peaks found.
thr_collapse_het_ratio = 0.1

```

## Build
To build from source.
```bash
git clone git@github.com:logsdon-lab/NucFreq.git && cd NucFreq
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
nucfreq -i test/HG00096_hifi_test.bam -b test/test.bed
```
```
haplotype2-0000133:3021508-8691473      3314093 3324276 COLLAPSE_VAR
haplotype2-0000133:3021508-8691473      4126277 4142368 MISJOIN
haplotype2-0000133:3021508-8691473      4566798 4683011 GAP
haplotype2-0000133:3021508-8691473      5466129 5496995 COLLAPSE
haplotype2-0000133:3021508-8691473      5737835 5747246 COLLAPSE_VAR
haplotype2-0000133:3021508-8691473      6067838 6072601 MISJOIN
haplotype2-0000133:3021508-8691473      6607947 6639102 COLLAPSE_VAR
haplotype2-0000133:3021508-8691473      7997560 8069465 MISJOIN
```

## Cite
- **Vollger MR**, Dishuck PC, Sorensen M, Welch AE, Dang V, Dougherty ML, et al. Long-read sequence and assembly of segmental duplications. Nat Methods. 2019;16: 88–94. doi:10.1038/s41592-018-0236-3
- **Mc Cartney AM**, Shafin K, Alonge M, Bzikadze AV, Formenti G, Fungtammasan A, et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. bioRxiv. 2021. p. 2021.07.02.450803. doi:10.1101/2021.07.02.450803
  * Citing `hetDetection.R`

## TODO
- Add false duplication detection.
- Publish on pypi.
- Refactor so cleaner.
- Colormap for `Misassembly`
