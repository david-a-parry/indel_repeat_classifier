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
(e.g. `indel_repeat_classifier --help`). Note that indel context is assessed
after simplifying and left-aligning variants. However the position, ref allele
and alt allele will be reported the same as provided in the VCF provided. You
may find it convenient to left-align variants prior to running
indel_repeat_classifier (e.g. using `bcftools norm`).

### Examples

Output repeat context information for indels in a given VCF using the
indel_repeat_classifier script:

    indel_repeat_classifier <in.vcf> <ref.fasta> -o output.csv -f 5

The output will be CSV format:

    chrom,pos,ref,alt,n_alleles,variant_length,cosmic_class,variant_type,repeat_type,repeat_unit,repeat_length,sequence
    1,25,CAAGG,C,2,4,4:Del:R:0,Del,No repeat,,0,agctgagagtcgtctcctcctcctcAAGGtcgtgcacagtctattgcacgtcgatgcgatgcgatgttg
    1,4,TGA,T,2,2,2:Del:R:1,Del,Perfect,GA,4,agctGAgagtcgtctcctcctcctca
    1,4,TGAGA,T,2,4,4:Del:M:1,Del,Imperfect,G,5,agctGAGAgtcgtctcctcctcctcaaggtcgtgcacagtctattgca
    1,13,TCTC,T,2,3,3:Del:R:3,Del,Perfect,CTC,12,agctgagagtcgtCTCctcctcctcaaggtcgtgcacagtctattg
    1,13,TCTCCTCCTCCTC,T,2,12,5:Del:R:0,Del,No repeat,,0,agctgagagtcgtCTCCTCCTCCTCaaggtcgtgcacagtctattgcacgtcgatgcgatgcgatgttgacagttagacacagtacacagtagagacagtag
    1,43,AT,A,2,1,1:Del:T:1,Del,Perfect,T,2,gcacagtctaTtgcacgtcga
    1,51,TCGATGCGATG,T,2,10,5:Del:M:5,Del,Imperfect,CGATG,15,agctgagagtcgtctcctcctcctcaaggtcgtgcacagtctattgcacgtCGATGCGATGcgatgttgacagttagacacagtacacagtagagacagtag
    1,51,TCGATG,T,2,5,5:Del:R:2,Del,Perfect,CGATG,15,agctgagagtcgtctcctcctcctcaaggtcgtgcacagtctattgcacgtCGATGcgatgcgatgttgacagttagacacagtacacagtagagacagtag
    1,27,A,ATC,2,2,2:Ins:R:0,Ins,No repeat,,0,agtcgtctcctcctcctcaaTCggtcgtgcacagtctattgc
    1,4,T,TGA,2,2,2:Ins:R:2,Ins,Perfect,GA,4,agctGAgagagtcgtctcctcctcct
    1,51,T,TCGATG,2,5,5:Ins:R:3,Ins,Perfect,CGATG,15,agctgagagtcgtctcctcctcctcaaggtcgtgcacagtctattgcacgtCGATGcgatgcgatgcgatgttgacagttagacacagtacacagtagagacagta
    1,26,A,AAGG,2,3,3:Ins:R:1,Ins,Perfect,AGG,3,agctgagagtcgtctcctcctcctcaAGGaggtcgtgcacagtctattgcacgtcgatg
    1,43,A,AT,2,1,1:Ins:T:2,Ins,Perfect,T,2,gcacagtctaTttgcacgtcg

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
[ID83](https://cancer.sanger.ac.uk/signatures/id/) signatures are also
provided.

By default, for variant data the Perfect/Imperfect classification refers to the
COSMIC ID83 style classification where, for example, a 4bp deletion of 'AGAG'
in a 6bp repeat of 'AGAGAG' is by default given the class of "Imperfect" repeat
as the deleted nucleotides are partially repeated. However, if you prefer to
classify this type of deletion as lying in a 'Perfect' repeat (because the
deletion is within a perfect 'AG' repeat) use the --collapse_alleles flag to
classify as a deletion at an "Perfect" repeat and to output the repeat unit as
'AG'. This will also classify deletions of entire repeat units as deletions of
a "Perfect" repeat as long as the deletion is limited to the repeat sequence.
Note that ID83 classifications will not be affected by this option.


## AUTHOR

Written by David A. Parry at the University of Edinburgh.

