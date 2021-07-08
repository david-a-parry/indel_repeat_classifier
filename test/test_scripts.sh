#!/bin/bash

set -euo pipefail
# quick test that scripts run as expected
bcftools concat test/test_data/variants/*vcf | \
    bin/indel_repeat_classifier - test/test_data/test.fasta

bcftools concat test/test_data/variants/*vcf | \
    bin/indel_repeat_classifier - test/test_data/test.fasta -f 5 -o test/deleteme.csv

bcftools concat test/test_data/variants/*vcf | \
    bin/indel_repeat_classifier - test/test_data/test.fasta -f 5 -c -o test/deleteme.csv.gz

bin/short_repeats_from_fasta test/test_data/test.fasta test/deleteme


echo $(date) Finished

echo "You may want to run 'rm test/deleteme*' now or after inspecting output files."
