#!/usr/bin/env python

"""Test pbtranscript.io.ReadStatRecord, ReadStatReader, ReadStatWriter."""
import unittest
import os.path as op
from pbtranscript.Utils import rmpath, mkdir
from pbtranscript.io.ReadStatIO import get_len_from_read_name, ReadStatRecord, ReadStatReader, ReadStatWriter
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR

_SIV_DIR_ = op.join(SIV_DATA_DIR, "test_ReadStatIO")
_DAT_DIR_ = op.join(DATA_DIR, "test_ReadStatIO")
_OUT_DIR_ = op.join(OUT_DIR, "test_ReadStatIO")

READSTAT_FN = op.join(_SIV_DIR_, "readstat.txt")

class TEST_ReadStatIO(unittest.TestCase):
    """Test classes pbtranscript.io.ReadStatRecord, ReadStatReader, ReadStatWriter."""
    def setUp(self):
        """Define input and output file."""
        rmpath(_OUT_DIR_)
        mkdir(_OUT_DIR_)

    def test_get_len_from_read_name(self):
        """Test get_len_from_read_name"""
        self.assertTrue(get_len_from_read_name("movie/124/100_300_CCS"), 200)
        self.assertTrue(get_len_from_read_name("movie/124/500_200_CCS"), 300)

    def test_ReadStatRecord(self):
        """test ReadSTatRecord."""
        record = ReadStatRecord(name="movie/123/123_0_CCS", is_fl=True,
                                stat="unique", pbid="PB.1.1")
        expected_str = "movie/123/123_0_CCS\t123\tY\tunique\tPB.1.1"
        self.assertEqual(str(record), expected_str)
        self.assertEqual(record.length, 123)

        expected_r = ReadStatRecord.fromString(expected_str)
        self.assertEqual(expected_r, record)

    def test_ReadStatReader(self):
        """"""
        reader = ReadStatReader(READSTAT_FN)
        records = [r for r in reader]
        self.assertEqual(len(records), 3)
        reader.close()
        expected_record_0 = "\t".join(["movie/100/123_0_CCS", "123", "Y", "unique", "PB.1.1"])
        expected_record_1 = "\t".join(["movie1/200/399_0_CCS", "399", "N", "unmapped", "NA"])
        expected_record_2 = "\t".join(["movie2/300/1_3991_CCS", "3990", "N", "ambiguous", "PB.2.1"])

        self.assertEqual(str(records[0]), expected_record_0)
        self.assertEqual(str(records[1]), expected_record_1)
        self.assertEqual(str(records[2]), expected_record_2)

        self.assertTrue(records[0].is_uniquely_mapped)
        self.assertTrue(records[1].is_unmapped)
        self.assertTrue(records[2].is_ambiguously_mapped)

        self.assertTrue(records[0].is_fl)
        self.assertFalse(records[1].is_fl)
        self.assertFalse(records[2].is_fl)

        out_fn = op.join(_OUT_DIR_, 'readstat.txt')

        writer = ReadStatWriter(out_fn, 'w')
        for r in records:
            writer.writeRecord(r)
        writer.close()
        records = [r for r in ReadStatReader(out_fn)]
        self.assertEqual(len(records), 3)

        writer = ReadStatWriter(out_fn, 'a')
        for r in records:
            writer.writeRecord(r)
        writer.close()
        records = [r for r in ReadStatReader(out_fn)]
        self.assertEqual(len(records), 3 * 2)
