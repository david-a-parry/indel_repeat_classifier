import re
from collections import namedtuple

nt_conversion = {'A': 'T', 'C': 'C', 'G': 'C', 'T': 'T'}

RepeatResult = namedtuple(
    "RepeatResult",
    'variant_type repeat_type repeat_unit repeat_length sequence')


def simplify_repeat(rpt):
    '''
    Represent a repeat unit in its simplest form.

    e.g. if regex found 'TTTTT' repeat unit, represent as 'T'
         if regex found 'TCTC' repeat unit, represent as 'TC'
    '''
    r_len = len(rpt)
    for i in range(1, int(r_len / 2) + 1):
        if r_len % i == 0:
            if rpt[:i] * int(r_len / i) == rpt:
                return rpt[:i]
    return rpt


def find_perfect_repeats_deletion(deletion, seq):
    '''
    Look for repeats of deletion at start of seq. Return length of repeat in bp
    '''
    match = re.match(r'(' + deletion + r')(\1)+', seq)
    if match:
        return match.end()
    return 0


def find_perfect_repeats_insertion(insertion, seq):
    '''
    Look for one or more occurences of insertion at start of seq.
    Return length of repeat in bp
    '''
    match = re.match(r"(" + insertion + r")+", seq)
    if match:
        return match.end()
    return 0


def find_microhomology(indel, seq, p):
    ''' Find microhomologies of indel at position p of seq. '''
    rpt_len, mh = 0, None
    i_len = len(indel)
    j = p + i_len
    for i in range(1, i_len):
        if seq[j:j + i] == indel[:i]:
            rpt_len = i
            mh = indel[:i]
        else:
            break
    indel = indel[::-1]
    seq = seq[::-1]
    p = len(seq) - p - i_len
    j = p + i_len
    for i in range(1, i_len):  # check reverse orientation
        if seq[j:j + i] == indel[:i]:
            if i > rpt_len:
                rpt_len = i
                mh = indel[:i]
        else:
            break
    if rpt_len:
        rpt_len += i_len
    return rpt_len, mh


def left_align_deletion(chrom, pos, ref, alt, fasta):
    ref_len = len(ref)
    alt_len = len(alt)
    if alt != ref[:alt_len]:
        return ref, alt, pos  # ignore complex indels
    seq = fasta[chrom]
    start = pos - 1
    end = start + ref_len
    while seq[start] == seq[end - 1] and start > 0:
        start -= 1
        end -= 1
        ref = seq[start:end]
        alt = ref[:alt_len]
    pos = start + 1
    return ref, alt, pos


def left_align_insertion(chrom, pos, ref, alt, fasta):
    ref_len = len(ref)
    alt_len = len(alt)
    if alt[:ref_len] != ref:
        return ref, alt, pos  # ignore complex indels
    seq = fasta[chrom]
    start = pos - 1
    end = start + ref_len
    var_len = alt_len - ref_len
    while alt[-1] == seq[start] and start > 0:
        start -= 1
        end -= 1
        alt = seq[start] + alt[:-1]
        ref = seq[start:end]
    pos = start + 1
    return ref, alt, pos


def left_align_variant(chrom, pos, ref, alt, fasta):
    ''' Left align a simplified variant to reference genome. '''
    ref_len = len(ref)
    alt_len = len(alt)
    if ref_len > alt_len and ref[:alt_len] == alt:
        return left_align_deletion(chrom, pos, ref, alt, fasta)
    elif ref_len < alt_len and alt[:ref_len] == ref:
        return left_align_insertion(chrom, pos, ref, alt, fasta)
    return ref, alt, pos


def simplify_variant(variant, fasta, allele=1):
    ''' Get the simplest representations of REF/ALT alleles and left-align. '''
    ref = variant.ref
    alt = variant.alleles[allele]
    pos = variant.pos
    while len(ref) > 1 and len(alt) > 1:
        if ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
        else:
            break
    while len(ref) > 1 and len(alt) > 1:
        if ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            pos += 1
        else:
            break
    return left_align_variant(variant.chrom, pos, ref, alt, fasta)


def repeats_from_variant(variant,
                         fasta,
                         allele=1,
                         min_flanks=10,
                         collapse_repeat_units=False):
    '''
    Find perfect repeats or microhomologies of indels. Assumes record comes
    from a normalized and left-aligned VCF.

    Returns a RepeatResult namedtuple with features 'variant_type',
    'repeat_type', 'repeat_unit', 'repeat_length' and 'sequence'. The
    repeat_length attribute gives either length of repeat or length of deletion
    plus microhomology in bp. The sequence attribute gives flanking sequence in
    lowercase and inserted/deleted bases in uppercase.

    Args:
        variant: pysam.VariantRecord

        fasta:  Fasta object from pyfaidx corresponding to reference sequence.

        allele: 1-based index of ALT allele to assess.

        min_flanks:
                Retrieve flanking sequence amounting to this number * variant
                length. For example, for a 2 bp deletion, setting this value to
                5 would result in retrieval and searching through 10 bp
                sequence either side of the variant. This value therefore caps
                the maximum repeat size that can be identified.

        collapse_repeat_units:
                Analyse repeats relative to the simplest repeating unit
                regardless of indel length. Default behaviour is to analyze
                repeats relative to the deleted/inserted nucleotides instead.

    '''
    var_type, rpt_type, rpt_unit, rpt_len, seq_ctxt = (None, "No repeat", None,
                                                       0, "")
    pos = variant.pos
    ref, alt, pos = simplify_variant(variant, fasta, allele)
    var_length = len(alt) - len(ref)
    flanks = abs(var_length) * min_flanks
    if var_length == 0:
        var_type = 'SNV' if len(ref) == 1 else 'MNV'
    elif (len(ref) != 1 and len(alt) != 1) or ref[0] != alt[0]:
        var_type = 'Complex'
    elif any(x not in 'ACGT' for x in ref + alt):
        var_type = 'Ambiguous nucleotide'
    if var_type is not None:
        return RepeatResult(var_type, rpt_type, rpt_unit, rpt_len, seq_ctxt)
    start = pos - 1
    stop = start + len(ref)
    l_flank = flanks - 1 if start - flanks > 0 else start
    r_flank = flanks if stop + flanks < len(fasta[variant.chrom]) \
        else len(fasta[variant.chrom]) - stop
    seq = fasta[variant.chrom][start - l_flank:stop + r_flank]
    seq_ctxt = seq[:l_flank].lower()
    if var_length < 0:
        var_type = 'Del'
        seq_ctxt += ref[0].lower() + ref[1:].upper() + \
            seq[l_flank + abs(var_length) + 1:].lower()
        if collapse_repeat_units:
            basic_rpt = simplify_repeat(ref[1:])
        else:
            basic_rpt = ref[1:]
        rpt_len = find_perfect_repeats_deletion(basic_rpt, seq[l_flank + 1:])
        if rpt_len:
            rpt_type = "Perfect"
            rpt_unit = basic_rpt
        else:
            rpt_len, rpt_unit = find_microhomology(ref[1:], seq, l_flank + 1)
            if rpt_len:
                rpt_type = "Imperfect"
    else:
        seq_ctxt += alt[0].lower() + alt[1:].upper() + \
            seq[l_flank + 1:].lower()
        var_type = 'Ins'
        if collapse_repeat_units:
            basic_rpt = simplify_repeat(alt[1:])
        else:
            basic_rpt = alt[1:]
        rpt_len = find_perfect_repeats_insertion(basic_rpt, seq[l_flank + 1:])
        if rpt_len:
            rpt_type = "Perfect"
            rpt_unit = basic_rpt
    return RepeatResult(var_type, rpt_type, rpt_unit, rpt_len, seq_ctxt)


def repeat_result_to_ID83(rpt_res, variant, fasta, allele):
    '''
    Convert RepeatResult (generated by repeat_from_variant) to COSMIC ID83
    classification.

    Args:
        rpt_res: RepeatResult namedtuple returned by repeats_from_variant (note
                 that repeats_from_variant should have been run with the
                 default collapse_repeat_units=False argument.

        variant: pysam.VariantRecord

        fasta:  Fasta object from pyfaidx corresponding to reference sequence.

        allele: 1-based index of ALT allele being assessed.

    '''
    rpt_size = 0
    var_len = abs(len(variant.ref) - len(variant.alleles[allele]))
    if rpt_res.repeat_type == 'Imperfect':
        return '{}:{}:M:{}'.format(min(var_len, 5), rpt_res.variant_type,
                                   min(rpt_res.repeat_length - var_len, 5))
    rpt_size = int(rpt_res.repeat_length / var_len)
    if rpt_size > 0 and rpt_res.variant_type == 'Del':
        rpt_size -= 1
    if var_len == 1:  # can only be perfect homopolymer repeat or no repeat
        ref, alt, pos = simplify_variant(variant, fasta, allele=allele)
        if len(variant.ref) > len(variant.alleles[allele]):
            nt = nt_conversion.get(ref[1], 'N')
        else:
            nt = nt_conversion.get(alt[1], 'N')
        return '1:{}:{}:{}'.format(rpt_res.variant_type, nt, min(rpt_size, 5))
    if rpt_size == 0 and rpt_res.repeat_type == 'Perfect':
        #  if entire deletion does not repeat COSMIC considers it microhomology
        #  even in case of e.g. 'GAGA' deletion in a perfect 'GAGAGA' repeat
        ref, alt, pos = simplify_variant(variant, fasta)
        start = pos - var_len
        end = pos + var_len * 2
        seq = fasta[variant.chrom][start:end]
        rpt_len, _ = find_microhomology(ref[1:], seq, var_len)
        if rpt_len:
            return '{}:Del:M:{}'.format(min(var_len, 5),
                                        min(rpt_len - var_len, 5))
    return '{}:{}:R:{}'.format(min(var_len, 5), rpt_res.variant_type,
                               min(rpt_size, 5))


def cosmic_ID83_classification(variant, fasta, allele=1):
    rpt_res = repeats_from_variant(variant, fasta, allele)
    return repeat_result_to_ID83(rpt_res, variant, fasta, allele)
