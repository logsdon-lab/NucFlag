```bash
zcat test/overlay/NA20847_rc-chr3_haplotype2-0000105\:89881870-96384969.bed.gz | awk -v OFS="\t" '{print $1, 0, $3}' | gzip > test/all_gap/cov.tsv.gz
```
