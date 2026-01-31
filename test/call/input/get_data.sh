#!/bin/bash

set -euox pipefail

WD=$(dirname "${0}")

infile_bam="/project/logsdon_shared/projects/Keith/NucFlagPaper/results/aligners/mm2/HG002_hifi.bam"
infile_fa="/project/logsdon_shared/data/reference/HG002/assembly/hg002v1.1.fasta"
bedfile="${WD}/HG002_chr1_MATERNAL.bed"

while read -r line; do
    region=$(printf "${line}" | awk '{ print $1":"$2"-"$3 }')
    region_fs=$(echo "${region}" | sed 's/:/_/g')
    outfile_bam="${WD}/HG002_${region_fs}.bam"
    outfile_fa="${WD}/HG002_${region_fs}.fa.gz"
    samtools view -bh "${infile_bam}" "${region}" -o "${outfile_bam}"
    samtools index "${outfile_bam}"

    samtools faidx "${infile_fa}" "${region}" | seqkit replace -p ":.+" | bgzip > "${outfile_fa}"
    samtools faidx "${outfile_fa}"
done < "${bedfile}"
