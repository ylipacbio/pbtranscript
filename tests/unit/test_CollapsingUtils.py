"""Test classes defined within pbtranscript.io.ContigSetReaderWrapper."""
import unittest
import os.path as op
from pbtranscript.Utils import rmpath
from pbtranscript.collapsing.CollapsingUtils import copy_sam_header, map_isoforms_and_sort
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR

_DIR_ = op.join(SIV_DATA_DIR, 'test_collapsing')

GMAP_INPUT_FASTA = op.join(_DIR_, 'gmap-input.fasta')
GMAP_INPUT_FASTQ = op.join(_DIR_, 'gmap-input.fastq')
GMAP_INPUT_FASTA_DS = op.join(_DIR_, 'gmap-input.fasta.contigset.xml')
GMAP_INPUT_FASTQ_DS = op.join(_DIR_, 'gmap-input.fastq.contigset.xml')
GMAP_SAM = op.join(SIV_DATA_DIR, 'test_SAMReader', 'gmap-output.sam')

GMAP_DB = op.join(SIV_DATA_DIR, 'gmap_db')
GMAP_NAME = 'SIRV'


class TEST_CollapsingUtils(unittest.TestCase):
    """Test functions of pbtranscript.collapsing.CollapsingUtils."""
    def setUp(self):
        """Define input and output file."""
        pass

    def test_copy_sam_header(self):
        """Test copy_sam_header"""
        out_fn = op.join(OUT_DIR, 'test_copy_sam_header.sam')
        copy_sam_header(GMAP_SAM, out_fn)
        with open(out_fn, 'r') as reader:
            lines = [r for r in reader]
            self.assertTrue(all([line.startswith('@') for line in lines]))
            self.assertEqual(len(lines), 9)

    def test_map_isoforms_and_sort(self):
        """Test map_isoforms_and_sort"""
        out_fn = op.join(OUT_DIR, 'test_map_isoforms_and_sort_fasta.sam')
        rmpath(out_fn)
        map_isoforms_and_sort(input_filename=GMAP_INPUT_FASTA,
                              sam_filename=out_fn,
                              gmap_db_dir=GMAP_DB,
                              gmap_db_name=GMAP_NAME,
                              gmap_nproc=10)
        self.assertTrue(op.exists(out_fn))

        out_fn = op.join(OUT_DIR, 'test_map_isoforms_and_sort_fastq.sam')
        rmpath(out_fn)
        map_isoforms_and_sort(input_filename=GMAP_INPUT_FASTQ,
                              sam_filename=out_fn,
                              gmap_db_dir=GMAP_DB,
                              gmap_db_name=GMAP_NAME,
                              gmap_nproc=10)
        self.assertTrue(op.exists(out_fn))

        out_fn = op.join(OUT_DIR, 'test_map_isoforms_and_sort_fasta_ds.sam')
        rmpath(out_fn)
        map_isoforms_and_sort(input_filename=GMAP_INPUT_FASTA_DS,
                              sam_filename=out_fn,
                              gmap_db_dir=GMAP_DB,
                              gmap_db_name=GMAP_NAME,
                              gmap_nproc=10)
        self.assertTrue(op.exists(out_fn))

        out_fn = op.join(OUT_DIR, 'test_map_isoforms_and_sort_fastq_ds.sam')
        rmpath(out_fn)
        map_isoforms_and_sort(input_filename=GMAP_INPUT_FASTQ_DS,
                              sam_filename=out_fn,
                              gmap_db_dir=GMAP_DB,
                              gmap_db_name=GMAP_NAME,
                              gmap_nproc=10)
        self.assertTrue(op.exists(out_fn))
