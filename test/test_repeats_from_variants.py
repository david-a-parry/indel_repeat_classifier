import sys
import os
import pysam
from pyfaidx import Fasta
from nose.tools import *
_t_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
_r_path = os.path.join(_t_path, os.pardir, 'repeat_metrics')
sys.path.insert(0, _r_path)
from repeats_from_variants import repeats_from_variant

dir_path = os.path.dirname(os.path.realpath(__file__))
var_path = os.path.join(dir_path, 'test_data', 'variants')
ref_fasta = os.path.join(dir_path, 'test_data', 'test.fasta')


def get_variants(path):
    records = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            records.append(rec)
    return records


def test_dels():
    vcf2expected = {"del1.vcf": ('Del', None, None, 0),
                    "del2.vcf": ('Del', 'Perfect', 'GA', 4),
                    "del3.vcf": ('Del', 'Perfect', 'GA', 4),
                    "del4.vcf": ('Del', 'Perfect', 'CTC', 12),
                    "del5.vcf": ('Del', 'Perfect', 'CTC', 12),
                    "del6.vcf": ('Del', 'Perfect', 'T', 2),
                    "del7.vcf": ('Del', 'Perfect', 'CGATG', 15),
                    }
    fasta = Fasta(ref_fasta, as_raw=True, sequence_always_upper=True)
    for vcf, expected in vcf2expected.items():
        records = get_variants(os.path.join(var_path, vcf))
        result = repeats_from_variant(records[0], fasta)
        assert_equal(result, expected)
    

def test_ins():
    vcf2expected = {"ins1.vcf": ('Ins', None, None, 0),
                    "ins2.vcf": ('Ins', 'Perfect', 'GA', 4),
                    "ins3.vcf": ('Ins', 'Perfect', 'CGATG', 15),
                    "ins4.vcf": ('Ins', 'Perfect', 'AGG', 3),
                    "ins5.vcf": ('Ins', 'Perfect', 'T', 2),
                    }
    fasta = Fasta(ref_fasta, as_raw=True, sequence_always_upper=True)
    for vcf, expected in vcf2expected.items():
        records = get_variants(os.path.join(var_path, vcf))
        result = repeats_from_variant(records[0], fasta)
        assert_equal(result, expected)
    

if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
