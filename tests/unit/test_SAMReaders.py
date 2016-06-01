#!/usr/bin/env python
"""
Test SAMReaders
"""

import unittest
import os.path as op

from pbtranscript.io.SAMReaders import GMAPSAMReader, SAMRecordBase, Interval, SAMflag

MNT_DATA = "/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data"
GMAP_SAM = op.join(MNT_DATA, "gmap-output.sam")

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


@unittest.skipUnless(op.isdir(MNT_DATA), "missing %s" % MNT_DATA)
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


if __name__ == "__main__":
    unittest.main()
