
# XXX verification for bug 30828 - runs ice_pbdagcon on a spurious cluster
# and checks that it discards the resulting all-N consensus sequence

import subprocess
import tempfile
import unittest
import os.path as op

from pbcore.io import FastaReader

CLUSTER_FA = """\
>m54007_151222_230824/47383194/334_64_CCS
CATTGAAGACGTCCACCTCAACGCTATGAACGTTAGTTGAGACAATGTTAAAGCAAACGACAACGTCATTGTGATCTACATACACAGTGGATGGTTAGCGTAAACATGGTGGAACGTACTTTGACTGCGCTGCAAGAAATGGTTGGGTCGATCGTAATGCTAGTCGTTACATCGGAACAAGCCAAAACAAAATCATTCGCTGGATTTAGACCTACTGCACGACGACGTCGACACAAGACATTCTTGAAAGGTAATTGACGTGGACGTTTC
>m54007_151222_230824/28640158/287_60_CCS
CAAACGACAACGTCATTGTGATCTACATACACAGTGGATGGTTAGGCGTAAACATGGTGGGAACGTACTTTGACTGCGCTGCAAGAAATGGGTTGGGTCGATCGTAATGCTAGTCGTTACATCGGAACAAGCCAAAAAACAAACATCATTCGCTGGATTTAGACTACTACTGCACGACCGACGTCGACACAAGACATTCTCTGAAAGGTAATTGACGTGGACGTTTC
>m54007_151222_230824/49611437/382_58_CCS
ACTGAACTACGGGTCAGCTTCCCCATTTGAAGTCATGTAGTGGTTGTCTACTTTTTCATTGAGACGTCCACCTCAACGCTATGAACGTTAGTTGAGACAATGTTAAAGCAAACGACAACGTCATTGTGATCTACATACACAGTGGATGGTTAGCGTAAACATGGTGGAACGTACTTTGACTGCGCTGCAAGAAATGGTGTGGGTCGATCGTAATGCTAGTCGTTACATCGGAACAAGCCAAAACAAAATCATTCGCTGGATTTAGACCTACTGCACGACGACGTCGACACAAGACATTCTTGAAAGGTAATTGACGTGGACGTT"""


class TestBadCluster(unittest.TestCase):

    def setUp(self):
        self.cluster_fa = tempfile.NamedTemporaryFile(suffix=".fasta").name
        with open(self.cluster_fa, "w") as fa_out:
            fa_out.write(CLUSTER_FA)

    def test_ice_pbdagcon_bad_cluster(self):
        out_fa = tempfile.NamedTemporaryFile(suffix=".fasta").name
        prefix = op.splitext(out_fa)[0]
        args = [
            "python", "-m", "pbtranscript.ice_pbdagcon",
            self.cluster_fa,
            prefix,
            "c5006"
        ]
        assert subprocess.call(args) == 0
        with FastaReader(out_fa) as fa_out:
            self.assertEqual(len([rec for rec in fa_out]), 0)


if __name__ == "__main__":
    unittest.main()
