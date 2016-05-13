
import unittest
import os.path as op

from pbtranscript.io.PbiBamIO import CCSInput

import test_setpath

DATA = test_setpath.DATA_DIR
print "DATA", DATA
DATA2 = "/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data/test_ccs_input/"


def _collect_sequence_stats(file_name):
    n_seqs = seq_len = file_len = 0
    with CCSInput(file_name) as f:
        file_len = len(f)  # test f.__len__()
        for rec in f:
            n_seqs += 1
            seq_len += len(rec.sequence)
    return n_seqs, seq_len, file_len


class TestCCSInput(unittest.TestCase):

    def test_fasta(self):
        file_name = op.join(DATA, "test.fasta")
        n_seqs, seq_len, file_len = _collect_sequence_stats(file_name)
        self.assertEqual(n_seqs, 3)
        self.assertEqual(seq_len, 6787)
        self.assertEqual(n_seqs, file_len)

    def test_contigset(self):
        file_name = op.join(DATA, "test_contigset.xml")
        n_seqs, seq_len, file_len = _collect_sequence_stats(file_name)
        self.assertEqual(n_seqs, 3)
        self.assertEqual(seq_len, 6787)
        self.assertEqual(n_seqs, file_len)

    @unittest.skipUnless(op.isdir(DATA2), "missing %s" % DATA2)
    def test_ccsset(self):
        file_name = op.join(DATA2, "tiny.ccs.xml")
        n_seqs, seq_len, file_len = _collect_sequence_stats(file_name)
        self.assertEqual(n_seqs, 14)
        self.assertEqual(seq_len, 14899)
        self.assertEqual(n_seqs, file_len)

if __name__ == "__main__":
    unittest.main()
