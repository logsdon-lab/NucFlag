
Format RM output:
```bash
 awk -v OFS="\t" 'NR >3 { print "NA20847_rc-chr3_haplotype2-0000105:89881871-96384969", $6, $7, $10, $11 }' /project/logsdon_share/projects/hgsvc3/CenMAP_verkko_run3/results/repeatmasker/NA20847_rc-chr3_haplotype2-0000105\:89881870-96384969_renamed.fa.out > test_cov/rm_renamed.bed
```

Correct coordinates.
```bash
awk -v OFS="\t" '{match($1, ":(.+)-", start); match($1, ".*-(.+)$", end); print $1, $2 + start[1], $3+start[1], $4, $5}' test_cov/rm_renamed.bed > test_cov/rm_renamed_adj.bed
```

Simplify annotations.
```bash
awk -v OFS="\t" '{
    split($5, rClass, "/" );
    new_rClass=rClass[1];
    if ($5 == "Satellite/centr" || $5 == "Satellite") {
        new_rClass=$4
    }
    switch (new_rClass) {
        case "SAR":
            new_rClass="HSat1A";
            break;
        case "HSAT":
            new_rClass="HSat1B";
            break;
        case "HSATII":
            new_rClass="HSat2";
            break;
        case "(CATTC)n":
            new_rClass="HSat2";
            break;
        case "(GAATG)n":
            new_rClass="HSat2";
            break;
        default:
            break;
    }
    print $1, $2, $3, new_rClass
}' test_cov/rm_renamed_adj.bed > test_cov/rm_renamed_adj_simple.bed
```

Reverse regions.
```bash
awk -v OFS="\t" '{new_start=96384969-$3+89881870;new_stop=96384969-$2+89881870; print $1, new_start, new_stop, $4, $5}' test/overlay/repeatmasker.bed | tac > test/overlay/repeatmasker_reort.bed
```
