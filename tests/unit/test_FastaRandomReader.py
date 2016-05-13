"""Test FastaSplitter"""

import unittest
import os.path as op
from pbcore.io import FastaReader, ContigSet
from pbtranscript.io.FastaRandomReader import FastaRandomReader, \
        MetaSubreadFastaReader
from pbtranscript.Utils import write_files_to_fofn
from test_setpath import DATA_DIR, OUT_DIR, STD_DIR
import hashlib

class TestFastaRandomReader(unittest.TestCase):
    """Class for testing FastaRandomReader."""

    def setUp(self):
        """Set up dataDir, outDir, stdoutDir"""
        self.outDir = OUT_DIR
        self.stdoutDir = STD_DIR
        self.inFa = op.join(DATA_DIR, "reads_of_insert.fasta")
        self.inXml = op.join(DATA_DIR, "reads_of_insert.contigset.xml")

        self.inFa2 = op.join(DATA_DIR, "primers.fasta")
        self.inXml2 = op.join(DATA_DIR, "test_ContigSetReaderWrapper.contigset.xml")
        self.multiFns = [self.inFa2, self.inXml2]

    def testSingleInput(self):
        """Test FastaRandomReader.keys() and __getitem__.
        Either a single FASTA as input or a single ContigSet xml as input.
        """
        reads = [r for r in FastaReader(self.inFa)]
        names = [r.name for r in reads]
        seqs = [r.sequence[:] for r in reads]

        self.assertEqual(len(reads), 22)

        xml_reads = [r for r in ContigSet(self.inXml)]
        xml_names = [r.name for r in xml_reads]
        xml_seqs  = [r.sequence for r in reads]

        # Reads in self.inFa and self.inXml should be identical.
        self.assertEqual(sorted(names), sorted(xml_names))
        self.assertEqual(sorted(seqs),  sorted(xml_seqs))

        # Random access input fasta file.
        frr = FastaRandomReader(self.inFa)
        self.assertEqual(sorted(list(frr.keys())), sorted(list(set(names))))

        for r in reads:
            self.assertEqual(frr[r.name].name, r.name)

        for r in reads:
            self.assertEqual(frr[r.name].sequence, r.sequence[:])

        # Random access input xml file
        frr = FastaRandomReader(self.inXml)
        self.assertEqual(sorted(list(frr.keys())), sorted(list(set(names))))

        for r in reads:
            self.assertEqual(frr[r.name].name, r.name)

        for r in reads:
            self.assertEqual(frr[r.name].sequence, r.sequence[:])

    def testMultipleInputs(self):
        """Test FastaRandomReader.keys() and __getitem__.
        Multiple input files.
        """
        reads = []
        for fn in self.multiFns:
            reads.extend([r for r in ContigSet(fn)])

        names = [r.name for r in reads]
        seqs = [r.sequence[:] for r in reads]

        self.assertEqual(len(reads), 29)

        # FastaRandomReader taking multiple files as input
        frr = FastaRandomReader(*self.multiFns)

        self.assertEqual(sorted(list(frr.keys())), sorted(list(set(names))))

        for r in reads:
            self.assertEqual(frr[r.name].name, r.name)

        for r in reads:
            self.assertEqual(frr[r.name].sequence, r.sequence[:])


class TestMetaSubreadFastaReader(unittest.TestCase):
    """Class for testing MetaSubreadFastaReader."""
    def setUp(self):
        """Set up testDir,  outDir, stdoutDir"""
        self.dataDir = DATA_DIR
        self.outDir = OUT_DIR

        self.fa1 = op.join(self.dataDir, "test_meta_subreads_fasta_reader1.fasta")
        self.fa2 = op.join(self.dataDir, "test_meta_subreads_fasta_reader2.fasta")
        self.fofn = op.join(self.outDir, "test_meta_subreads_fasta_reader.fofn")

    def testAll(self):
        """Test FastaRandomReader.keys() and __getitem__."""
        write_files_to_fofn([self.fa1, self.fa2], self.fofn)
        reader = MetaSubreadFastaReader([self.fa1, self.fa2])
        subread_1 = "m130812_185809_42141_c100533960310000001823079711101380_s1_p0/59/0_5071"
        subread_2 = "m130812_random_random_s1_p0/440/13280_16126"
        zmw_3 = "m130812_185809_42141_c100533960310000001823079711101380_s1_p0/70"
        zmw_4 = "m130812_random_random_s1_p0/249"
        r1 = reader[subread_1][0]
        self.assertEqual(r1.name, subread_1)
        self.assertEqual(hashlib.md5(r1.sequence).hexdigest(), "8128261dd851ae285d029618739559e9")

        r2 = reader[subread_2][0]
        self.assertEqual(r2.name, subread_2)
        self.assertEqual(hashlib.md5(r2.sequence).hexdigest(), "451e5798a7f21cce80da27a03a8cb2c7")

        r3, r4 = reader[zmw_3]
        self.assertEqual(r3.name, "m130812_185809_42141_c100533960310000001823079711101380_s1_p0/70/0_5538")
        self.assertEqual(r4.name, "m130812_185809_42141_c100533960310000001823079711101380_s1_p0/70/5587_5982")
        self.assertEqual(hashlib.md5(r3.sequence).hexdigest(), "4db2f6e35c83dd279a8f71a51ac50445")
        self.assertEqual(hashlib.md5(r4.sequence).hexdigest(), "1c1d080e9362a73ea2074f9a62fbd45e")

        r5 = reader[zmw_4][0]
        self.assertEqual(r5.name, "m130812_random_random_s1_p0/249/0_1339")
        self.assertEqual(hashlib.md5(r5.sequence).hexdigest(), "b20d3723a136aedc2f96f6f498ad3da0")

