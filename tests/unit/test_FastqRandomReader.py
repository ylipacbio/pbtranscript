"""Test FastqRandomReader"""

import unittest
import os.path as op
from pbcore.io import FastqReader
from pbtranscript.io.FastqRandomReader import FastqRandomReader
from test_setpath import DATA_DIR


class TestFastqRandomReader(unittest.TestCase):
    """Class for testing FastqRandomReader."""
    def setUp(self):
        """Set up dataDir."""
        self.dataDir = DATA_DIR
        self.inFq = op.join(self.dataDir, "test_fastq_random_reader.fastq")

    def testAll(self):
        """Test FastqRandomReader.keys() and __getitem__."""
        reads = [r for r in FastqReader(self.inFq)]
        names = [r.name for r in reads]

        frr = FastqRandomReader(self.inFq)
        self.assertTrue(set(frr.keys()) == set(names))

        self.assertTrue(False not in
                [frr[r.name].name == r.name for r in reads])

        self.assertTrue(False not in
                [frr[r.name].sequence == r.sequence for r in reads])

        self.assertTrue(False not in
                [frr[r.name].quality.all() == r.quality.all() for r in reads])

        self.assertTrue(False not in
                [frr[r.name].qualityString == r.qualityString for r in reads])
