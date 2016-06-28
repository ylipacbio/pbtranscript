#!/usr/bin/env python

"""Test pbtranscript.io.AbundanceRecord, AbundanceReader, AbundanceWriter."""
import unittest
import os.path as op
from pbtranscript.io import AbundanceRecord, AbundanceReader, AbundanceWriter
import filecmp
from test_setpath import DATA_DIR, OUT_DIR

ABUNDANCE_FN = op.join(DATA_DIR, "test_Abundance.txt")

class TEST_AbundanceIO(unittest.TestCase):
    """Test classes pbtranscript.io.AbundanceRecord, AbundanceReader, AbundanceWriter."""
    def setUp(self):
        """Define input and output file."""
        pass

    def test_AbundanceRecord(self):
        """test AbundanceRecord."""
        record = AbundanceRecord(pbid="PB.1.1", count_fl=100, count_nfl=200, count_nfl_amb=300.00,
                                 norm_fl=0.2000, norm_nfl=0.3000, norm_nfl_amb=0.4000)
        expected_str = "PB.1.1\t100\t200\t300.00\t2.0000e-01\t3.0000e-01\t4.0000e-01"
        self.assertEqual(str(record), expected_str)

        new_record = AbundanceRecord.fromString(expected_str)
        self.assertEqual(str(new_record), expected_str)

    def test_AbundanceReader_Writer(self):
        """test AbundanceReader and AbundanceWriter"""
        reader = AbundanceReader(ABUNDANCE_FN)
        records = [r for r in reader]
        self.assertEqual(len(records), 3)
        reader.close()

        expected_total_fl = 4712
        expected_total_nfl = 9879
        expected_total_nfl_amb = 10064

        self.assertEqual(reader.total_fl, expected_total_fl)
        self.assertEqual(reader.total_nfl, expected_total_nfl)
        self.assertEqual(reader.total_nfl_amb, expected_total_nfl_amb)

        out_fn = op.join(OUT_DIR, "test_Abundance.txt")
        writer = AbundanceWriter(out_fn, reader.comments)
        for r in records:
            writer.writeRecord(r)
        writer.close()
        self.assertTrue(filecmp.cmp(out_fn, ABUNDANCE_FN))

        out_fn = op.join(OUT_DIR, "test_Abundance2.txt")
        writer = AbundanceWriter(out_fn, None, expected_total_fl,
                                 expected_total_nfl, expected_total_nfl_amb)
        for r in records:
            writer.writeRecord(r)
        writer.close()
        self.assertTrue(filecmp.cmp(out_fn, ABUNDANCE_FN))
