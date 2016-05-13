#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbcore.io import FastaReader
import pbcommand.testkit.core
from pbtranscript.Utils import execute, mknewdir
from pbtranscript.separate_flnc import SizeBin, SizeBins, \
        SeparateFLNCByPrimer, SeparateFLNCBySize
from test_setpath import DATA_DIR, OUT_DIR, STD_DIR, SIV_DATA_DIR

PBI_DATA = "/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data/bam"
FLNC_DATASET = op.join(PBI_DATA, "isoseq_flnc.contigset.xml")
FLNC_FASTA = op.join(PBI_DATA, "isoseq_flnc.fasta")


def test_SizeBin():
    """Test SizeBin"""
    b1 = SizeBin(1,2)
    b2 = SizeBin(2,3)
    b3 = SizeBin(3,4)

    print "b1=%s" % str(b1)
    print "b2=%s" % str(b2)
    print "%s" % (b1 < b2)
    #assert (b1 < b2)

    assert str(b1) == "1to2kb"
    assert b1.contains(1500) is True
    assert b1.contains(500) is False
    assert b1.contains(2500) is False

def test_SizeBins():
    """Test class SizeBins"""
    bs = SizeBins([0, 2, 3, 6])
    assert bs[0] == SizeBin(0, 2)
    assert bs[1] == SizeBin(2, 3)
    assert bs[2] == SizeBin(3, 6)

    assert bs.which_bin_contains(3900) == SizeBin(3, 6)
    assert bs.which_bin_contains(SizeBin(3,4)) == SizeBin(3, 6)


class TestSeparateFLNCByPrimer(unittest.TestCase):
    """Test SeparateFLNCByPrimer"""
    def setUp(self):
        """Initialize."""
        pass

    def test_run_fasta_in(self):
        """Test run()."""
        expected_sorted_keys = [1]
        out_dir=op.join(OUT_DIR, 'separate_flnc_by_primer_fasta_input')
        has_contigset_output=False

        with SeparateFLNCByPrimer(flnc_filename=FLNC_FASTA,
                                  root_dir=out_dir) as obj:
            obj.run()

        self.assertEqual(obj.sorted_keys, expected_sorted_keys)

        for index, key in enumerate(obj.sorted_keys):
            with FastaReader(obj.out_fasta_files[index]) as reader:
                self.assertTrue(all([obj._get_primer_id(r) == key for r in reader]))

        for xml_fn in obj.out_contigset_files:
            print xml_fn
            self.assertEqual(op.exists(xml_fn), has_contigset_output)

    def test_run_xml_in(self):
        """Test run(), contigset.xml as input."""
        expected_sorted_keys = [1]
        out_dir = op.join(OUT_DIR, 'separate_flnc_by_primer_xml_input')
        has_contigset_output = True

        with SeparateFLNCByPrimer(flnc_filename=FLNC_DATASET,
                                  root_dir=out_dir) as obj:
            obj.run()

        self.assertEqual(obj.sorted_keys, expected_sorted_keys)

        for index, key in enumerate(obj.sorted_keys):
            with FastaReader(obj.out_fasta_files[index]) as reader:
                self.assertTrue(all([obj._get_primer_id(r) == key for r in reader]))

        for xml_fn in obj.out_contigset_files:
            print xml_fn
            self.assertEqual(op.exists(xml_fn), has_contigset_output)



class TestSeparateFLNCBySize(unittest.TestCase):
    """Test SeparateFLNCBySize"""
    def setUp(self):
        """Initialize."""
        pass

    def _test_bin_manual(self, bin_manual, expected_bin_manual):
        """Test SeparateFLNCBySize setting bin manually."""
        out_dir=op.join(OUT_DIR, 'separate_flnc_by_size_bin_manual')
        mknewdir(out_dir)
        with SeparateFLNCBySize(flnc_filename=FLNC_FASTA,
                                bin_manual=bin_manual,
                                root_dir=out_dir) as obj:
             obj.run()

        self.assertEqual(obj.sorted_keys, expected_bin_manual)

        for index, key in enumerate(obj.sorted_keys):
            with FastaReader(obj.out_fasta_files[index]) as reader:
                self.assertTrue(all([key[0].contains(len(r.sequence)) for r in reader]))

    def test_bin_manual(self):
        """Test run()."""
        bin_manual = []
        print 'testing %s' % bin_manual
        expected_bin_manual = [(SizeBin(3, 4), 0), (SizeBin(4,5), 0)] # [(SizeBin, part_idx), ... ]
        self._test_bin_manual(bin_manual=bin_manual, expected_bin_manual=expected_bin_manual)

        bin_manual = [3, 10]
        print 'testing %s' % bin_manual
        expected_bin_manual = [(SizeBin(3, 10), 0)]
        self._test_bin_manual(bin_manual=bin_manual, expected_bin_manual=expected_bin_manual)

        bin_manual = [0, 1]
        print 'testing %s' % bin_manual
        expected_bin_manual = [(SizeBin(1, 5), 0)]
        self._test_bin_manual(bin_manual=bin_manual, expected_bin_manual=expected_bin_manual)

        bin_manual = [2, 3, 4, 8]
        print 'testing %s' % bin_manual
        expected_bin_manual = [(SizeBin(3, 4), 0), (SizeBin(4, 8), 0)]
        self._test_bin_manual(bin_manual=bin_manual, expected_bin_manual=expected_bin_manual)


def test_end_to_end():
    """Call separate_flnc.py from command line, end to end must exit gracefully."""
    cmd = "separate_flnc.py %s %s --bin_by_primer" % \
          (FLNC_FASTA, op.join(OUT_DIR, "separate_flnc_by_primer_fasta_input_e2e"))
    execute(cmd)

    cmd = "separate_flnc.py %s %s --bin_by_primer" % \
          (FLNC_DATASET, op.join(OUT_DIR, "separate_flnc_by_primer_xml_input_e2e"))
    execute(cmd)

    cmd = "separate_flnc.py %s %s --bin_size_kb 1" % \
          (FLNC_FASTA, op.join(OUT_DIR, "separate_flnc_by_size_fasta_input_e2e"))
    execute(cmd)

    cmd = "separate_flnc.py %s %s --bin_size_kb 1" % \
          (FLNC_DATASET, op.join(OUT_DIR, "separate_flnc_by_size_xml_input_e2e"))
    execute(cmd)

    cmd = "separate_flnc.py %s %s --bin_manual '[0,3,4,6]'" % \
          (FLNC_FASTA, op.join(OUT_DIR, "separate_flnc_by_size_fasta_input_manual_e2e"))
    execute(cmd)

    cmd = "separate_flnc.py %s %s --bin_manual '[0,3,4,6]'" % \
          (FLNC_DATASET, op.join(OUT_DIR, "separate_flnc_by_size_xml_input_manual_e2e"))
    execute(cmd)
