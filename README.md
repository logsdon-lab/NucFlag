# NucFreq Misassembly Classifier
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

## Build
To build from source.
```bash
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
make test
```

Or try the test example.
```bash
nucfreq -i test/HG00096_hifi_test.bam -b test/test/bed
```

## Cite
- **Vollger MR**, Dishuck PC, Sorensen M, Welch AE, Dang V, Dougherty ML, et al. Long-read sequence and assembly of segmental duplications. Nat Methods. 2019;16: 88–94. doi:10.1038/s41592-018-0236-3
- **Mc Cartney AM**, Shafin K, Alonge M, Bzikadze AV, Formenti G, Fungtammasan A, et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. bioRxiv. 2021. p. 2021.07.02.450803. doi:10.1101/2021.07.02.450803
  * Citing `hetDetection.R`

## TODO
- Add false duplication detection.
- Publish on pypi.
