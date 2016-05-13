"""Test FastaSplitter"""

import unittest
import os.path as op
from pbcore.io import FastaReader
from pbtranscript.io.FastaSplitter import FastaSplitter

class TestFastaSplitter(unittest.TestCase):
    """Class for testing FastaSplitter."""
    def setUp(self):
        """Set up testDir, dataDir, out_dir, stdout_dir"""
        self.rootDir = op.dirname(op.dirname(op.abspath(__file__)))
        self.testDir = op.join(self.rootDir, "")
        self.dataDir = op.join(self.testDir, "data")
        self.out_dir = op.join(self.testDir, "out")
        self.stdout_dir = op.join(self.testDir, "stdout")
        self.input_fasta = op.join(self.dataDir, "reads_of_insert.fasta")

    def testSplit(self):
        """Test FastaSplitter.split()."""
        fs = FastaSplitter(self.input_fasta, 2, self.out_dir,
            "testFastaSplitter_split_")
        fs.split()
        splittedReads = []
        for of in fs.out_fns:
            self.assertTrue(op.exists(of))
            with FastaReader(of) as reader:
                splittedReads.extend([(r.name, r.sequence) for r in reader])
        fs.rmOutFNs()

        reads = []
        with FastaReader(self.input_fasta) as reader:
            reads.extend([(r.name, r.sequence) for r in reader])
        self.assertTrue(len(reads) == 22)
        self.assertTrue(splittedReads == reads)

    def testRmOutFNs(self):
        """Test FastaSplitter.rmOutFNs()."""
        fs = FastaSplitter(self.input_fasta, 2, self.out_dir,
            "testFastaSplitter_rmOutFNs_")
        fs.split()
        fs.rmOutFNs()
        self.assertTrue(len(fs.out_fns) == 0)




