#!/usr/bin/env python
"""
Test SAMReaders
"""

import unittest
import os.path as op

from pbtranscript.io.GffIO import GffRecordBase, CollapseGffRecord, \
        CollapseGffReader, CollapseGffWriter, GmapRecord
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR


GFF_FN = op.join(SIV_DATA_DIR, 'test_GffIO', "collapse.gff")


class TestGffRecordBase(unittest.TestCase):

    def test_all(self):
        """Test GffRecordBase"""
        record = GffRecordBase(seqid="seqid", source=None,
                               feature="feature", start=0, end=100,
                               score=None, strand="+", frame=None,
                               attributes="gene_id \"PB.1\"; transcript_id \"PB.1.1\";")
        expected_record_str = "seqid\t.\tfeature\t0\t100\t.\t+\t.\tgene_id \"PB.1\"; transcript_id \"PB.1.1\";"
        self.assertEqual(str(record), expected_record_str)

        other = GffRecordBase.fromString(expected_record_str)
        self.assertEqual(str(other), expected_record_str)
        self.assertEqual(record, other)

@unittest.skipUnless(op.isdir(SIV_DATA_DIR), "Missing %s" % SIV_DATA_DIR)
class TestCollapseGffRecord(unittest.TestCase):
    def test_all(self):
        """Test CollapseGffRecord"""
        record = CollapseGffRecord(seqid="seqid", feature="exon",
                                   start=0, end=100,
                                   strand="+", gene_id="PB.1",
                                   transcript_id="PB.1.1")
        expected_record_str = "seqid\tPacBio\texon\t0\t100\t.\t+\t.\tgene_id \"PB.1\"; transcript_id \"PB.1.1\";"
        self.assertEqual(str(record), expected_record_str)

        other = CollapseGffRecord.fromString(expected_record_str)
        self.assertEqual(str(other), expected_record_str)
        self.assertEqual(record, other)

        self.assertTrue(record.is_exon)
        self.assertFalse(record.is_transcript)

    def test_reader_writer(self):
        """Test CollapseGffReader."""
        expected_num_records = 58
        expected_num_exons_in_records = [3, 6, 3, 6, 4, 5, 3, 11, 11, 2, 5, 3, 3, 5, 3, 5, 2, 8, 3, 3, 3, 3, 2, 2, 4, 5, 16, 16, 2, 17, 2, 17, 6, 18, 17, 17, 4, 1, 9, 10, 10, 10, 8, 9, 4, 4, 4, 5, 5, 6, 1, 3, 5, 6, 5, 3, 5, 6]
        expected_total_num_exons = 422 - expected_num_records
        with CollapseGffReader(GFF_FN) as reader:
            records = [r for r in reader]
            self.assertEqual(len(records), expected_num_records)
            self.assertEqual([len(r.ref_exons) for r in records], expected_num_exons_in_records)
            self.assertEqual(sum(expected_num_exons_in_records), expected_total_num_exons)

        with CollapseGffReader(GFF_FN) as reader:
            records = []
            r = None
            while True:
                try:
                    records.append(reader.next())
                except StopIteration:
                    break
            self.assertEqual(len(records), expected_num_records)
            self.assertEqual([len(r.ref_exons) for r in records], expected_num_exons_in_records)
            self.assertEqual(sum(expected_num_exons_in_records), expected_total_num_exons)

    def test_writer(self):
        """
        Test CollapseGffWriter.

        Read -> Write -> Read yields identical record.
        """
        out_fn = op.join(OUT_DIR, "test_GffIO.collapse.gff")
        expected_records = [r for r in CollapseGffReader(GFF_FN)]
        with CollapseGffWriter(out_fn) as writer:
            for r in expected_records:
                writer.writeRecord(r)

        records = [r for r in CollapseGffReader(out_fn)]
        expected_r0 = expected_records[0]
        r0 = records[0]
        self.assertEqual(records, expected_records)
