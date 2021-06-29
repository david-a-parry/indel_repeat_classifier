import sys
import os
import pysam
from pyfaidx import Fasta
_t_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
_r_path = os.path.join(_t_path, os.pardir, 'indel_repeat_classifier')
sys.path.insert(0, _r_path)
from repeats_from_variants import repeats_from_variant
from repeats_from_variants import cosmic_ID83_classification

dir_path = os.path.dirname(os.path.realpath(__file__))
var_path = os.path.join(dir_path, 'test_data', 'variants')
fa_path = os.path.join(dir_path, 'test_data', 'refs')
ref_fasta = os.path.join(dir_path, 'test_data', 'test.fasta')


def get_variants(path):
    records = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            records.append(rec)
    return records


def test_dels():
    ''' Identify deletions in perfect repeats '''
    vcf2expected = {
        'del1.vcf':
        ('Del', None, None, 0,
         'agctgagagtcgtctcctcc' + 'tcctcAAGGtcgtgcacagtctattgcacgtcg'),
        'del2.vcf': ('Del', 'Perfect', 'GA', 4, 'agctGAgagtcgtctcct'),
        'del3.vcf':
        ('Del', 'Perfect', 'GA', 4, 'agctGAGAgtcgtctcctcctcctcaaggtcg'),
        'del4.vcf': ('Del', 'Perfect', 'CTC', 12,
                     'agctgagagtcgtCTCctcctcctcaaggtcgtg'),
        'del5.vcf':
        ('Del', 'Perfect', 'CTC', 12,
         'agctgagagtcgtCTCCTCCTCCTCaaggtcgtgcacagtct' +
         'attgcacgtcgatgcgatgcgatgttgacagttagacacagt' + 'acacagtagagac'),
        'del6.vcf': ('Del', 'Perfect', 'T', 2, 'agtctaTtgcacg'),
        'del7.vcf':
        ('Del', 'Perfect', 'CGATG', 15,
         'agctgagagtcgtctcctcctcctcaaggtcgtgcacagtct' +
         'attgcacgtCGATGCGATGcgatgttgacagttagacacagt' + 'acacagtagagacagtag'),
    }
    fasta = Fasta(ref_fasta, as_raw=True, sequence_always_upper=True)
    for vcf, expected in vcf2expected.items():
        records = get_variants(os.path.join(var_path, vcf))
        result = repeats_from_variant(records[0], fasta, min_flanks=6)
        assert result == expected


def test_ins():
    ''' Identify insertions in perfect repeats '''
    vcf2expected = {
        "ins1.vcf": ('Ins', None, None, 0,
                     'agtcgtctcctcctcctca' + 'aTCggtcgtgcacagtctattgc'),
        "ins2.vcf": ('Ins', 'Perfect', 'GA', 4, 'agctGAgagagtcgtctcctcctcct'),
        "ins3.vcf": ('Ins', 'Perfect', 'CGATG', 15,
                     'agctgagagtcgtctcctcctcctcaaggtcgtgcacagtct' +
                     'attgcacgtCGATGcgatgcgatgcgatgttgacagttagac' +
                     'acagtacacagtagagacagta'),
        "ins4.vcf":
        ('Ins', 'Perfect', 'AGG', 3,
         'agctgagagtcgtctcctcctcctcaAGGaggtcgtgcacag' + 'tctattgcacgtcgatg'),
        "ins5.vcf": ('Ins', 'Perfect', 'T', 2, 'gcacagtctaTttgcacgtcg'),
    }
    fasta = Fasta(ref_fasta, as_raw=True, sequence_always_upper=True)
    for vcf, expected in vcf2expected.items():
        records = get_variants(os.path.join(var_path, vcf))
        result = repeats_from_variant(records[0], fasta)
        assert result == expected


def test_microhomology():
    ''' Identify deletions with microhomology '''
    size2del = {
        2: ('TC', 'cccccccccccccccaccaaTCtagcggcccccccccccccc'),
        3: ('TTC', 'c' * 25 + 'acccaTTCtagcgg' + 'c' * 24,
            'c' * 25 + 'acccaTTCttagcgg' + 'c' * 23),
        4: ('TATC', 'c' * 35 + 'acccaTATCttagcgg' + 'c' * 33,
            'c' * 35 + 'acccaTATCtaagcgg' + 'c' * 33,
            'c' * 35 + 'acccaTATCtatagcgg' + 'c' * 32),
        5: ('TAGTC', 'c' * 45 + 'acccaTAGTCttagcgg' + 'c' * 43,
            'c' * 45 + 'acccaTAGTCtaagcgg' + 'c' * 43,
            'c' * 45 + 'acccaTAGTCtagagcgg' + 'c' * 42,
            'c' * 45 + 'acccaTAGTCtagtagcgg' + 'c' * 41),
        7: ('TAGCCTC', 'c' * 50 + 'acccaTAGCCTCtagcctagcgg' + 'c' * 47)
    }
    for i in range(2, 6):
        for j in range(i - 1, i):
            base = "mh_{}_{}_f".format(i, j)
            mh_fa = os.path.join(fa_path, base + '.fasta')
            records = get_variants(os.path.join(var_path, base + '.vcf'))
            fasta = Fasta(mh_fa, as_raw=True, sequence_always_upper=True)
            mh_seq = size2del[i][0][:j]
            seq_context = size2del[i][j]
            expected = ('Del', 'Imperfect', mh_seq, i + j, seq_context)
            result = repeats_from_variant(records[0], fasta)
            assert result == expected
    for base, mh_seq, j, sc in zip(
        ["mh_2_1_r", "mh_7_6_f"], ['C', 'TAGCCT'], [3, 13], [
            'c' * 16 + 'acaaCTcaagcgg' + 'c' * 13,
            'c' * 50 + 'acccaTAGCCTCtagcctagcgg' + 'c' * 47
        ]):
        mh_fa = os.path.join(fa_path, base + '.fasta')
        records = get_variants(os.path.join(var_path, base + '.vcf'))
        fasta = Fasta(mh_fa, as_raw=True, sequence_always_upper=True)
        expected = ('Del', 'Imperfect', mh_seq, j, sc)
        result = repeats_from_variant(records[0], fasta)
        assert result == expected


def test_del_cosmic_classification():
    ''' Test COSMIC ID83 classification of deletions'''
    vcf2expected = {
        "del1.vcf": '4:Del:R:0',
        "del2.vcf": '2:Del:R:1',
        "del3.vcf": '4:Del:M:1',
        "del4.vcf": '3:Del:R:3',
        "del5.vcf": '5:Del:R:0',
        "del6.vcf": '1:Del:T:1',
        "del7.vcf": '5:Del:M:5',
    }
    fasta = Fasta(ref_fasta, as_raw=True, sequence_always_upper=True)
    for vcf, expected in vcf2expected.items():
        records = get_variants(os.path.join(var_path, vcf))
        result = cosmic_ID83_classification(records[0], fasta)
        expected = vcf2expected[vcf]
        assert result == expected


def test_mh_cosmic_classification():
    ''' Test COSMIC ID83 classification of deletions with microhomology'''
    for i in range(2, 6):
        for j in range(i - 1, i):
            base = "mh_{}_{}_f".format(i, j)
            mh_fa = os.path.join(fa_path, base + '.fasta')
            records = get_variants(os.path.join(var_path, base + '.vcf'))
            fasta = Fasta(mh_fa, as_raw=True, sequence_always_upper=True)
            expected = '{}:Del:M:{}'.format(i, j)
            result = cosmic_ID83_classification(records[0], fasta)
            assert result == expected


def test_ins_cosmic_classification():
    ''' Test COSMIC ID83 classification of insertions'''
    vcf2expected = {
        "ins1.vcf": '2:Ins:R:0',
        "ins2.vcf": '2:Ins:R:2',
        "ins3.vcf": '5:Ins:R:3',
        "ins4.vcf": '3:Ins:R:1',
        "ins5.vcf": '1:Ins:T:2'
    }
    fasta = Fasta(ref_fasta, as_raw=True, sequence_always_upper=True)
    for vcf, expected in vcf2expected.items():
        records = get_variants(os.path.join(var_path, vcf))
        result = cosmic_ID83_classification(records[0], fasta)
        assert result == expected


if __name__ == '__main__':
    import nose2
    nose2.main()
