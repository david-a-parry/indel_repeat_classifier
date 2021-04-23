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
    if [[ "$var" =~ mh_[0-9]_[0-9]_ ]]
    then
        ref_fasta=${source_dir}/refs/$(basename $var )
    else
        ref_fasta=$ref
    fi
    bowtie2 -D 20 -R 10  -N 1 -L 5  --rdg 1,1 --rfg 1,1  -f -x $ref_fasta -U $var \
        --rg "ID:Sample1" --rg "SM:Sample1" | \
        bcftools mpileup --no-version --fasta-ref $ref_fasta - | \
        bcftools call --no-version -P 1 -m -v -O u - | \
        bcftools norm --no-version -f $ref_fasta > ${var_dir}/$(basename $var .fasta).vcf
done

cp -v ${source_dir}/complex.vcf ${var_dir}/

echo $(date) Done
echo $?
