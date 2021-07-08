import csv
import gzip
import os
import sys
import tempfile

_t_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
_r_path = os.path.join(_t_path, os.pardir, 'indel_repeat_classifier')
sys.path.insert(0, _r_path)
from repeats_from_sequence import find_repeats_in_fasta

dir_path = os.path.dirname(os.path.realpath(__file__))
csv_path = os.path.join(dir_path, 'test_data', 'csvs')
ref_fasta = os.path.join(dir_path, 'test_data', 'test.fasta')


def test_repeats_in_fasta():
    "Identify short repeats in FASTA sequences"
    for i in range(1, 6):
        exp_csv = os.path.join(csv_path, "{}bp_repeats.csv".format(i))
        seq_name = "{}bp_rpts".format(i)
        _, tmp_path = tempfile.mkstemp(suffix='.csv.gz')
        find_repeats_in_fasta(ref_fasta,
                              tmp_path,
                              rpt_unit_length=i,
                              flanks=5,
                              seqnames=[seq_name],
                              quiet=True)
        with open(exp_csv, 'rt') as exp_fh, gzip.open(tmp_path,
                                                      'rt') as res_fh:
            exp_reader = csv.DictReader(exp_fh)
            res_reader = csv.DictReader(res_fh)
            for exp in exp_reader:
                res = next(res_reader)
                assert res == exp


if __name__ == '__main__':
    import nose2
    nose2.main()
