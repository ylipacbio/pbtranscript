
"""
Test for appropriate handling of different input types in the low-level IO
code used by ICE.
"""

import unittest
import os.path as op

from pbcore.io import FastaReader

from pbtranscript.ice.ProbModel import ProbFromQV
from pbtranscript.io.BasQV import basQVcacher

MNT_DATA = "/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data"
CCS_BAM = MNT_DATA + "/movies/rat_bax1/Analysis_Results/m131018_081703_42161_c100585152550000001823088404281404_s1_p0.1.ccs.bam"
CCS_XML = op.join(MNT_DATA, "sa3", "ccsbam.consensusreadset.xml")
CCS_FOFN = op.join(MNT_DATA, "ccsbam.fofn")
NFL_FASTA = op.join(MNT_DATA, "nfl.fasta")

READ_ID = "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/43/0_568_CCS"

def _get_read_ids():
    with FastaReader(NFL_FASTA) as f:
        for read in f:
            yield read.name.split()[0]


@unittest.skipUnless(op.isdir(MNT_DATA), "missing %s" % MNT_DATA)
class TestBasQVCacher(unittest.TestCase):

    def test_fofn(self):
        pm = ProbFromQV(op.join(MNT_DATA, CCS_FOFN),
            fasta_filename=NFL_FASTA)
        qvs = []
        for read_id in _get_read_ids():
            qvs.append(pm.get(read_id, "InsertionQV"))
        dqv = pm.get(READ_ID, "DeletionQV")
        self.assertEqual("%.5f" % dqv[0], "0.01995")
        #print dqv[100]
        self.assertEqual(len(qvs), 251)

    def test_xml(self):
        pm = ProbFromQV(op.join(MNT_DATA, CCS_XML),
            fasta_filename=NFL_FASTA)
        qvs = []
        for read_id in _get_read_ids():
            qvs.append(pm.get(read_id, "InsertionQV"))
        dqv = pm.get(READ_ID, "DeletionQV")
        self.assertEqual("%.5f" % dqv[0], "0.01995")
        #print dqv[100]
        self.assertEqual(len(qvs), 251)

    def test_bam(self):
        qver = basQVcacher()
        qver.add_bash5(CCS_BAM)
        seqids = [ rid for rid in _get_read_ids() ]
        qver.precache(seqids)
        qvs = []
        for read_id in seqids:
            qvs.append(qver.get(read_id, "InsertionQV"))
        dqv = qver.get(READ_ID, "DeletionQV")
        self.assertEqual("%.5f" % dqv[0], "0.01995")
        #print dqv[100]
        self.assertEqual(len(qvs), 251)


if __name__ == "__main__":
    unittest.main()
