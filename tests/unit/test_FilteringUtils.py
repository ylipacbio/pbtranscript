"""Test pbtranscript.filtering.FilteringUtils."""
import unittest
import os.path as op
import cPickle
import filecmp
import numpy as np

from pbcore.io import FastqReader
from pbtranscript.io import CollapseGffReader, AbundanceReader, GroupReader
from pbtranscript.Utils import rmpath, mkdir
from pbtranscript.filtering.FilteringUtils import good_isoform_ids_by_count, \
    good_isoform_ids_by_removing_subsets, filter_by_count, filter_out_subsets

from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR, SIV_STD_DIR

GROUP_FN = op.join(SIV_DATA_DIR, "test_filtering", "in.group.txt")
ABUNDANCE_FN = op.join(SIV_DATA_DIR, "test_filtering", "in.abundance.txt")
GFF_FN = op.join(SIV_DATA_DIR, "test_filtering", "in.gff")
REP_FN = op.join(SIV_DATA_DIR, "test_filtering", "in.rep.fastq")

_OUT_DIR_ = op.join(OUT_DIR, "test_filtering")
rmpath(_OUT_DIR_)
mkdir(_OUT_DIR_)


class TEST_FilteringUtils(unittest.TestCase):
    """Test functions of pbtranscript.filtering.FilteringUtils."""
    def setUp(self):
        """Define input and output file."""
        self.expected_good = ['PB.2.5', 'PB.5.1', 'PB.7.1', 'PB.10.2', 'PB.10.42', 'PB.12.1']
        self.expected_diff = ['PB.10.42', 'PB.10.36', 'PB.10.35']

    def test_good_isoform_ids_by_count(self):
        """Test good_isoform_ids_by_count"""
        good = good_isoform_ids_by_count(in_group_filename=GROUP_FN,
                                         in_abundance_filename=ABUNDANCE_FN,
                                         min_count=20)
        self.assertEqual(good, self.expected_good)

    def test_good_isoform_ids_by_removing_subsets(self):
        """Test good_isoform_ids_by_removing_subsets"""
        all = [r.seqid for r in CollapseGffReader(GFF_FN)]
        good = good_isoform_ids_by_removing_subsets(in_gff_filename=GFF_FN,
                max_fuzzy_junction=5)
        diff = list(set(all) - set(good))
        self.assertEqual(diff, self.expected_diff)

    def test_filter_by_count(self):
        """Test filter_by_count"""
        out_abundance_fn = op.join(_OUT_DIR_, "filter_by_count.abundance.txt")
        out_gff_fn = op.join(_OUT_DIR_, "filter_by_count.gff")
        out_rep_fn = op.join(_OUT_DIR_, "filter_by_count.rep.fastq")
        filter_by_count(in_group_filename=GROUP_FN, in_abundance_filename=ABUNDANCE_FN,
                        in_gff_filename=GFF_FN, in_rep_filename=REP_FN,
                        out_abundance_filename=out_abundance_fn,
                        out_gff_filename=out_gff_fn,
                        out_rep_filename=out_rep_fn,
                        min_count=20)

        out_abundance_ids = [r.pbid for r in AbundanceReader(out_abundance_fn)]
        self.assertEqual(out_abundance_ids, self.expected_good)

        out_gff_ids = [r.seqid for r in CollapseGffReader(out_gff_fn)]
        self.assertEqual(out_gff_ids, self.expected_good)

        out_rep_ids = [r.name.split('|')[0] for r in FastqReader(out_rep_fn)]
        self.assertEqual(out_rep_ids, self.expected_good)

    def test_filter_out_subsets(self):
        """Test filter_out_subsets"""
        out_abundance_fn = op.join(_OUT_DIR_, "filter_out_subsets.abundance.txt")
        out_gff_fn = op.join(_OUT_DIR_, "filter_out_subsets.gff")
        out_rep_fn = op.join(_OUT_DIR_, "filter_out_subsets.rep.fastq")
        filter_out_subsets(in_abundance_filename=ABUNDANCE_FN,
                           in_gff_filename=GFF_FN, in_rep_filename=REP_FN,
                           out_abundance_filename=out_abundance_fn,
                           out_gff_filename=out_gff_fn, out_rep_filename=out_rep_fn,
                           max_fuzzy_junction=5)

        all = [r.seqid for r in CollapseGffReader(GFF_FN)]
        expected_good = set(all)-set(self.expected_diff)
        out_abundance_ids = [r.pbid for r in AbundanceReader(out_abundance_fn)]
        self.assertEqual(set(out_abundance_ids), expected_good)

        out_gff_ids = [r.seqid for r in CollapseGffReader(out_gff_fn)]
        self.assertEqual(set(out_gff_ids), expected_good)

        out_rep_ids = [r.name.split('|')[0] for r in FastqReader(out_rep_fn)]
        self.assertEqual(set(out_rep_ids), expected_good)
