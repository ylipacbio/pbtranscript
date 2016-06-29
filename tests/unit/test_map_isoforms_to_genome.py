"""Test pbtranscript.collapsing.Branch."""
import unittest
import os.path as op
import cPickle
import filecmp
import numpy as np
from pbtranscript.Utils import rmpath, mkdir
from pbtranscript.tasks.map_isoforms_to_genome import gmap_db_and_name_from_ds
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR, SIV_STD_DIR


READS_DS = op.join(SIV_DATA_DIR, 'test_collapsing', 'gmap-input.fastq.contigset.xml')
GMAP_DS = op.join(SIV_DATA_DIR, "gmap-referenceset-root-dir/SIRV/gmapreferenceset.xml")
_OUT_DIR_ = op.join(OUT_DIR, "test_map_isoforms_to_genome")


class TEST_map_isoforms_to_genome(unittest.TestCase):
    """Test functions of pbtranscript.tasks.map_isoforms_to_genome."""
    def setUp(self):
        """Define input and output file."""
        rmpath(_OUT_DIR_)
        mkdir(_OUT_DIR_)

    def test_gmap_db_and_name_from_ds(self):
        """Test map_isoforms_to_genome.gmap_db_and_name_from_ds"""
        gmap_db, gmap_name = gmap_db_and_name_from_ds(GMAP_DS)
        self.assertEqual(gmap_db, op.join(SIV_DATA_DIR, "gmap-referenceset-root-dir", "SIRV"))
        self.assertEqual(gmap_name, "gmap_db")

