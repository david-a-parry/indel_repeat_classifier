# Indel Repeat Classifier

Scripts and python modules for assessing short simple repeats in either the
context of variants (insertions/deletions) or a reference sequence.

## INSTALLATION

Python3 is required. `pysam` and `pyfaidx` modules are also required and can be
installed using `pip`.

Due to the requirement for the
[pysam](https://pysam.readthedocs.io/en/latest/api.html) module only linux or
Mac OS X systems are supported.

To install directly from github:

    pip3 install git+git://github.com/david-a-parry/indel_repeat_classifier.git --user

Alternatively, you may first clone this repository:

    git clone https://github.com/david-a-parry/indel_repeat_classifier

From the newly created indel_repeat_classifier directory you may install either
by running the setup.py script as follows:

    python3 setup.py install --user

or by using pip, if installed:

    pip3 install . --user

## USAGE

For all available options for the provided scripts run with the `--help` flag
(e.g. `indel_repeat_classifier --help`). Please note that any VCF provided to
the indel_repeat_classifier script must first be normalized and left-aligned
(e.g. using `bcftools norm`).

### Examples

Output repeat context information for indels in a given VCF using the
indel_repeat_classifier script:

    indel_repeat_classifier <in.vcf> <ref.fasta> -o output.csv -f 5

The output will be CSV format:

    chrom,pos,ref,alt,n_alleles,variant_length,cosmic_class,variant_type,repeat_type,repeat_unit,repeat_length,sequence,cosmic_repeat_type
    1,25,CAAGG,C,2,4,4:Del:R:0,Del,No repeat,,0,agagtcgtctcctcctcctcAAGGtcgtgcacagtctattgcac,No repeat
    1,4,TGA,T,2,2,2:Del:R:1,Del,Perfect,GA,4,agctGAgagtcgtctc,Perfect
    1,4,TGAGA,T,2,4,4:Del:M:1,Del,Perfect,GA,4,agctGAGAgtcgtctcctcctcctcaag,Imperfect
    1,13,TCTC,T,2,3,3:Del:R:3,Del,Perfect,CTC,12,agctgagagtcgtCTCctcctcctcaaggtc,Perfect
    1,13,TCTCCTCCTCCTC,T,2,12,5:Del:R:0,Del,Perfect,CTC,12,agctgagagtcgtCTCCTCCTCCTCaaggtcgtgcacagtctattgcacgtcgatgcgatgcgatgttgacagttagacacagta,No repeat

Output all repeats between of 2-4 bp in a given fasta file:

    short_repeats_from_fasta <fasta> <output_prefix> -r 2 3 4

Output files will be named <output_prefix>\_<repeat_unit_size>.csv.gz

Alternatively, to use as a module for classifying variants:

    #!/usr/bin/env python3
    import pysam
    from pyfaidx import Fasta
    from indel_repeat_classifier.repeats_from_variants import repeats_from_variant

    ref_fasta = Fasta("/path/to/genome.fa")
    with pysam.VariantFile("path/to/vcf") as variants:
        for record in variants:
            for i in range(len(record.alleles)):
                rpt_res = repeats_from_variant(variant=record,
                                               fasta=ref_fasta,
                                               allele=i)
                # do something with repeat result
                ...

## OUTPUT

Repeat types are classified as "Perfect" for complete repeats (e.g. "AGAG" is a
2bp perfect repeat of "AG") while partial repeats (sometimes referred to as
repeats with microhomology) are classified as "Imperfect" (e.g. "AGAT" is an
imperfect repeat of "AG"). Sequence context is given with the deleted bases (or
repeat unit if using `short_repeats_from_fasta`) in uppercase and flanking bases
in lowercase. For variant data the cosmic classifications relating to the
[ID83](https://cancer.sanger.ac.uk/signatures/id/) signatures are also provided.

By default, for variant data the Perfect/Imperfect classification refers to
whether the repeat overlapping the deletion is a Perfect/Imperfect repeat , but
this behaviour can be modified using the `--classify_allele` flag where the
Perfect/Imperfect classification refers to whether the insertion/deletion allele
itself is repeated. For example, a 4bp deletion of 'AGAG' in a 6bp repeat of
'AGAGAG' is by default given the class of "Perfect" repeat as it lies within a
perfect 'AG' repeat. However, if you prefer the SigProfiler interpretation
where this would be considered microhomology rather than a deletion of a perfect
perfect repeat (e.g. "4:Del:M:2") use this flag to classify as a deletion at an
"Imperfect" repeat.
