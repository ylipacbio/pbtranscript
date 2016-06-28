#!/usr/bin/env python
"""
Test SAMReaders
"""

import unittest
import os.path as op

from pbtranscript.Utils import rmpath
from pbtranscript.io import ContigSetReaderWrapper
from pbtranscript.io.SAMReaders import GMAPSAMReader, SAMRecordBase, Interval, SAMflag, iter_gmap_sam
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR

GMAP_SAM = op.join(SIV_DATA_DIR, 'test_SAMReader', "gmap-output.sam")
SORTED_GMAP_SAM = op.join(SIV_DATA_DIR, 'test_SAMReader', "sorted-gmap-output.sam")
READS_DS = op.join(SIV_DATA_DIR, 'test_SAMReader', 'gmap-input.fastq')


def _get_sam_groups(ignored_ids_writer=None):
    """Returns grouped sam records read from SORTED_GMAP_SAM and READS_DS."""
    query_len_dict = ContigSetReaderWrapper.name_to_len_dict(READS_DS)
    groups = [g for g in iter_gmap_sam(sam_filename=SORTED_GMAP_SAM,
                                       query_len_dict=query_len_dict,
                                       min_aln_coverage=0.99, min_aln_identity=0.85,
                                       ignored_ids_writer=ignored_ids_writer)]
    return groups


def construct_SAMRecord(qID, sID, qStart, qEnd, sStart, sEnd, qLen, sLen,
                        cigar, flag, segments,
                        mismatches, ins, de, mat_or_sub,
                        identity, qCov, sCov):
    """Return a SAMRecordBase object."""
    obj = SAMRecordBase()
    obj.qID, obj.sID = qID, sID
    obj.qStart, obj.qEnd = qStart, qEnd
    obj.sStart, obj.sEnd = sStart, sEnd
    obj.qLen, obj.sLen = qLen, sLen
    obj.cigar, obj.flag, obj.segments = cigar, flag, segments
    obj.num_nonmatches = mismatches
    obj.num_ins, obj.num_del, obj.num_mat_or_sub = ins, de, mat_or_sub
    obj.identity, obj.qCoverage, obj.sCoverage = identity, qCov, sCov
    return obj


@unittest.skipUnless(op.isdir(SIV_DATA_DIR), "missing %s" % SIV_DATA_DIR)
class TestSAMReaders(unittest.TestCase):

    def test_header(self):
        """Test GMAPReader.header"""
        expected_ref_len_dict = dict({'SIRV5': 14606, 'SIRV4': 16122,
                                      'SIRV7': 148957, 'SIRV6': 12837,
                                      'SIRV1': 12643, 'SIRV3': 10943,
                                      'SIRV2': 6911})
        reader = GMAPSAMReader(GMAP_SAM)
        self.assertEqual(reader.header.referenceLengthsDict,
                         expected_ref_len_dict)


    def test_sam_records(self):
        reads = [r for r in GMAPSAMReader(GMAP_SAM)]
        self.assertEqual(len(reads), 984)

        expected_r0 = construct_SAMRecord(
                "i0_HQ_sampleb92221|c24/f8p0/933", "SIRV5",
                0, 933, 12670, 13603,
                None, 14606,
                "933M",
                SAMflag(False, '+', 0),
                [Interval(start=12670, end=13603)],
                0, 0, 0, 933,
                1.0, None, 0.0638778584144)
        expected_r973 = construct_SAMRecord(
                "i2_HQ_sampleb92221|c4709/f3p0/2081", "SIRV5",
                0, 2081, 1009, 10987,
                None, 14606,
                "50M1I20M1D69M838N46M86N37M114N45M983N106M79N160M1737N70M93N83M485N58M158N125M206N64M104N131M187N163M374N81M108N511M73N131M2273N84M2I9M1D35M",
                SAMflag(False, '+', 0),
                [Interval(start=1009, end=1149), Interval(start=1987, end=2033),
                 Interval(start=2119, end=2156), Interval(start=2270, end=2315),
                 Interval(start=3298, end=3404), Interval(start=3483, end=3643),
                 Interval(start=5380, end=5450), Interval(start=5543, end=5626),
                 Interval(start=6111, end=6169), Interval(start=6327, end=6452),
                 Interval(start=6658, end=6722), Interval(start=6826, end=6957),
                 Interval(start=7144, end=7307), Interval(start=7681, end=7762),
                 Interval(start=7870, end=8381), Interval(start=8454, end=8585),
                 Interval(start=10858, end=10987)],
                5, 3, 2, 2078,
                0.997599615939, None, 0.68314391346)
        expected_r981 = construct_SAMRecord(
                "i2_HQ_sampleb92221|c4719/f2p0/2242", "SIRV7",
                0, 2242, 1003, 114915,
                None, 148957,
                "1672M318N118M984N84M546N85M38218N49M71603N180M1D54M",
                SAMflag(False, '-', 0),
                [Interval(start=1003, end=2675), Interval(start=2993, end=3111),
                 Interval(start=4095, end=4179), Interval(start=4725, end=4810),
                 Interval(start=43028, end=43077), Interval(start=114680, end=114915)],
                1, 0, 1, 2242,
                0.999554168524, None, 0.764730761226)

        self.assertEqual(reads[0], expected_r0)
        self.assertEqual(reads[973], expected_r973)
        self.assertEqual(reads[981], expected_r981)


    def test_iter_gmap_sam(self):
        """
        test iter_gmap_sam, which takes a sorted gmap sam file as input, and
        iterates over a group of overlapping sam records (which supposed to belong
        to the same isoform family.)
        """
        ignored_ids_txt = op.join(OUT_DIR, 'iter_gmap_sam.ignored.txt')
        rmpath(ignored_ids_txt)
        ignored_ids_writer = open(ignored_ids_txt, 'w')
        groups = _get_sam_groups(ignored_ids_writer)
        ignored_ids_writer.close()

        self.assertTrue(op.exists(ignored_ids_txt))
        ignored_ids = [line.split(' ')[0] for line in open(ignored_ids_txt, 'r')]
        self.assertEqual(len(ignored_ids), 108)

        self.assertEqual(len(groups), 9)
        expected_plus_lens = [10, 2, 129, 31, 0, 0, 348, 141, 0]
        self.assertEqual([len(g["+"]) for g in groups], expected_plus_lens)

        expected_minus_lens = [77, 36, 11, 0, 6, 9, 2, 2, 72]
        self.assertEqual([len(g["-"]) for g in groups], expected_minus_lens)

        self.assertTrue(all([r.sID == 'SIRV1' for r in groups[0]["+"]]))
        self.assertTrue(all([r.sID == 'SIRV2' for r in groups[1]["+"]]))
        self.assertTrue(all([r.sID == 'SIRV3' for r in groups[2]["+"]]))
        self.assertTrue(all([r.sID == 'SIRV4' for r in groups[3]["+"]]))
        self.assertTrue(all([r.sID == 'SIRV4' for r in groups[4]["-"]]))
        self.assertTrue(all([r.sID == 'SIRV4' for r in groups[5]["-"]]))
        self.assertTrue(all([r.sID == 'SIRV5' for r in groups[6]["+"]]))
        self.assertTrue(all([r.sID == 'SIRV6' for r in groups[7]["+"]]))
        self.assertTrue(all([r.sID == 'SIRV7' for r in groups[8]["-"]]))

        expected_g0_plus_sStart = [10710, 10712, 10712, 10712, 10712, 10712, 10712, 10713, 10713, 10715]
        expected_g0_plus_sEnd = [11641, 11641, 11638, 11640, 11641, 11641, 11638, 11641, 11640, 11641]
        self.assertTrue(expected_g0_plus_sStart, [r.sStart for r in groups[0]["+"]])
        self.assertTrue(expected_g0_plus_sEnd, [r.sEnd for r in groups[0]["+"]])

        expected_g4_minus_sStart = [3640, 3640, 3642, 3642, 3642, 3644]
        expected_g4_minus_sEnd = [5157, 5157, 5157, 5157, 3829, 5157]
        self.assertTrue(expected_g0_plus_sStart, [r.sStart for r in groups[4]["-"]])
        self.assertTrue(expected_g0_plus_sEnd, [r.sEnd for r in groups[4]["-"]])


if __name__ == "__main__":
    unittest.main()
