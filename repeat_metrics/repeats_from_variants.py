import re

nt_conversion = {'A': 'T',
                 'C': 'C',
                 'G': 'G',
                 'T': 'T'}

def simplify_repeat(rpt):
    '''
    Represent a repeat unit in its simplest form.

    e.g. if regex found 'TTTTT' repeat unit, represent as 'T'
         if regex found 'TCTC' repeat unit, represent as 'TC'
    '''
    r_len = len(rpt)
    for i in range(1, int(r_len/2) + 1):
        if r_len % i == 0:
            if rpt[:i] * int(r_len/i) == rpt:
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
        if seq[j:j+i] == indel[:i]:
            rpt_len = i
            mh = indel[:i]
        else:
            break
    indel = indel[::-1]
    seq = seq[::-1]
    p = len(seq) - p - i_len
    j = p + i_len
    for i in range(1, i_len):  # check reverse orientation
        if seq[j:j+i] == indel[:i]:
            if i > rpt_len:
                rpt_len = i
                mh = indel[:i]
        else:
            break
    return rpt_len, mh


def simplify_variant(variant, allele=1):
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
    return ref, alt, pos


def repeats_from_variant(variant, fasta, allele=1, min_flanks=20):
    '''
    Find perfect repeats or microhomologies of indels. Assumes record comes
    from a normalized and left-aligned VCF.

    Return variant type, type of repeat, repeat_unit and either length
    of repeat in units or length of microhomology in bp.

    Args:
        variant: pysam.VariantRecord

        fasta:  Fasta object from pyfaidx corresponding to reference sequence.

        allele: 1-based index of ALT allele to assess.

    '''
    var_type, rpt_type, rpt_unit, rpt_len = None, None, None, 0
    pos = variant.pos
    ref, alt, pos = simplify_variant(variant)
    var_length = len(alt) - len(ref)
    if var_length == 0:
        var_type = 'SNV' if len(ref) == 1 else 'MNV'
    if (len(ref) != 1 and len(alt) != 1) or ref[0] != alt[0]:
        var_type = 'Complex'
    if var_type is not None:
        return var_type, rpt_type, rpt_len
    start = pos - 1
    stop = start + len(ref)
    flanks = min_flanks if min_flanks > abs(var_length) \
        else abs(var_length) + min_flanks
    l_flank = flanks if start - flanks > 0 else start
    r_flank = flanks if stop + flanks < len(fasta[variant.chrom]) \
        else len(fasta[variant.chrom]) - stop
    seq = fasta[variant.chrom][start - l_flank:stop + r_flank]
    if var_length < 0:
        var_type = 'Del'
        basic_rpt = simplify_repeat(ref[1:])
        rpt_len = find_perfect_repeats_deletion(basic_rpt, seq[l_flank + 1:])
        if rpt_len:
            rpt_type = "Perfect"
            rpt_unit = basic_rpt
        else:
            rpt_len, rpt_unit = find_microhomology(ref[1:], seq, l_flank + 1)
            if rpt_len:
                rpt_type = "Imperfect"
    else:
        var_type = 'Ins'
        basic_rpt = simplify_repeat(alt[1:])
        rpt_len = find_perfect_repeats_insertion(basic_rpt, seq[l_flank + 1:])
        if rpt_len:
            rpt_type = "Perfect"
            rpt_unit = basic_rpt
    return var_type, rpt_type, rpt_unit, rpt_len


def cosmic_ID83_classification(variant, fasta, allele=1):
    var_type, rpt_type, rpt_unit, rpt_len = repeats_from_variant(variant,
                                                                 fasta,
                                                                 allele)
    rpt_size = 0
    var_len = abs(len(variant.ref) - len(variant.alleles[allele]))
    if rpt_type == 'Imperfect':
        return '{}:{}:M:{}'.format(min(var_len, 5), var_type, rpt_len)
    rpt_size = int(rpt_len/var_len)
    if rpt_size > 0 and var_type == 'Del':
        rpt_size -= 1
    rpt_size = rpt_size if rpt_size < 5 else 5
    if var_len == 1:  # can only be perfect homopolymer repeat or no repeat
        ref, alt, pos = simplify_variant(variant)
        if var_type == 'Del':
            nt = nt_conversion[ref[1]]
        else:
            nt = nt_conversion[alt[1]]
        return ('1:{}:{}:{}'.format(var_type, nt, rpt_size))
    return '{}:{}:R:{}'.format(min(var_len, 5), var_type, min(rpt_size, 5))
