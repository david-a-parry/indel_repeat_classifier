#!/bin/bash

source_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ref=${source_dir}/test.fasta
var_dir=${source_dir}/variants
set -exou pipefail

mkdir -p $var_dir

pyfasta split --header="${var_dir}/%(seqid)s.fasta" ${source_dir}/indels.fasta 
rename -v s/\ .+/.fasta/ ${var_dir}/*.fasta

bowtie2-build $ref $ref

for var in ${var_dir}/*fasta
do 
    echo $(date) Processing $var
    bowtie2 -D 20 -R 3 -N 1 -L 5  --rdg 1,1 --rfg 1,1  -f -x $ref -U $var \
        --rg "ID:Sample1" --rg "SM:Sample1" | \
        bcftools mpileup --fasta-ref $ref - | \
        bcftools call -P 1 -m -v -O u - | \
        bcftools norm -f $ref > ${var_dir}/$(basename $var .fasta).vcf
done

cp -v ${source_dir}/complex.vcf ${var_dir}/

echo $(date) Done
echo $?
