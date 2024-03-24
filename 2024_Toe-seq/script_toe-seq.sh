#!/bin/bash

# Script for obtaining group_by_seqs.txt file from PE reads
# arg1 = R1.fq
# arg2 = R2.fq

set -eo pipefail

cutadapt -g GACTCGAGTCGACATCGA -a CAGCCGTATGGATTTAACACG -G CGTGTTAAATCCATACGGCTG -A TCGATGTCGACTCGAGTC --report=minimal --match-read-wildcards --cores=13 --trim-n --times=10 --discard-untrimmed --pair-filter=both --overlap=7 --minimum-length=10 -o trimmed1_R1.fq -p trimmed1_R2.fq $1 $2
cutadapt -g CGTGTTAAATCCATACGGCTG -a TCGATGTCGACTCGAGTC -G GACTCGAGTCGACATCGA -A CAGCCGTATGGATTTAACACG --report=minimal --match-read-wildcards --cores=13 --trim-n --times=10 --discard-untrimmed --pair-filter=both --overlap=7 --minimum-length=10 -o trimmed2_R1.fq -p trimmed2_R2.fq $1 $2

cat trimmed1_R1.fq trimmed2_R1.fq >trimmed_R1.fq
cat trimmed1_R2.fq trimmed2_R2.fq >trimmed_R2.fq
rm trimmed1_R1.fq trimmed2_R1.fq trimmed1_R2.fq trimmed2_R2.fq

$usearch -fastq_mergepairs trimmed_R1.fq -reverse trimmed_R2.fq -fastqout merged.fq -fastq_minovlen 10 -minhsp 10 -fastq_maxdiffs 5 -fastq_pctid 90 -fastq_trunctail 25 -fastq_minlen 10 -sample ${1%_R1*}
$usearch -fastq_filter merged.fq -fastaout merged.fa -fastq_maxee_rate 0.005
rm trimmed_R1.fq trimmed_R2.fq merged.fq

cutadapt -g TTTTTTTTTTTT...CAGAACCCCACC -g TTTTTTTTTTTT...AAGCAGGAAGCC -g TTTTTTTTTTTT...GGAAGGCGATAT -g TTTTTTTTTTTT...CGCAAGCCGCCG --report=minimal --overlap 8 --discard-untrimmed --minimum-length 21 --cores=13 --action=none -o with_adapt.fa merged.fa
rm merged.fa

cutadapt -g GCATTAATACGACTCACTATAGGG -e 0.15 --report=minimal --overlap 13 --minimum-length 21 --cores=13 --untrimmed-output wo_T7.fa -o with_T7.fa with_adapt.fa
rm with_adapt.fa with_T7.fa

python3 scripts/trim_fasta_by_ref.py -f wo_T7.fa -b -d scripts/db.txt -e 1 -n scripts/db_const.txt --split CCCACC >with_barcodes1.fa
python3 scripts/trim_fasta_by_ref.py -f wo_T7.fa -b -d scripts/db.txt -e 1 -n scripts/db_const.txt --split GAAGCC >with_barcodes2.fa

cat with_barcodes1.fa with_barcodes2.fa >with_barcodes.fa
rm wo_T7.fa with_barcodes1.fa with_barcodes2.fa

python3 scripts/group_seq_by_seq.py -f with_barcodes.fa -d scripts/db.txt -c scripts/db_const.txt --max 71 >group_by_seqs.txt
rm with_barcodes.fa

