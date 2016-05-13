"""Test pbtranscript.SubsetExtractor."""

import unittest
import os.path as op
from pbtranscript.SubsetExtractor import SubsetRules, \
    ReadsSubsetExtractor
from pbtranscript.io.ReadAnnotation import ReadAnnotation
from pbcore.io import FastaReader
import filecmp
from test_setpath import DATA_DIR, OUT_DIR, STD_DIR

class Test_SubsetExtractor(unittest.TestCase):
    """Test SubsetExtractor."""
    def setUp(self):
        """Set up test data."""
        self.dataDir = DATA_DIR
        self.outDir = OUT_DIR
        self.stdDir = STD_DIR

    def test_satisfy(self):
        """Test function satisfy()."""
        inFN = op.join(self.dataDir, "test_subset.fasta")
        reads = []
        with FastaReader(inFN) as reader:
            reads = [x for x in reader]

        rules = SubsetRules(1, 1) # Full-length, non-chimeric
        obj = ReadsSubsetExtractor("in", "out", rules, False, True)
        #inFN, outFN, rules, ignore_polyA, printReadLengthOnly):

        ans = [ReadAnnotation.fromString(r.name) for r in reads]
        res = [obj.satisfy(an, rules) for an in ans]
        expected = [1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1]
        self.assertTrue(res == expected)

    def test_run(self):
        """Test function run()."""
        inFN = op.join(self.dataDir, "test_subset.fasta")
        outFN = op.join(self.outDir, "test_subset_unit.fasta")
        stdoutFN = op.join(self.stdDir, "test_subset_unit.fasta")

        rules = SubsetRules(1, 1) # Full-length, non-chimeric
        obj = ReadsSubsetExtractor(inFN, outFN, rules, False, True)
        obj.run()
        self.assertTrue(filecmp.cmp(outFN, stdoutFN))





