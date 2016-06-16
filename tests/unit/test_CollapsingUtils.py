"""Test classes defined within pbtranscript.collapsing.CollaspingUtils."""
import unittest
import os.path as op
from pbtranscript.Utils import rmpath, mkdir
from pbtranscript.io import CollapseGffReader, CollapseGffRecord
from pbtranscript.collapsing.CollapsingUtils import copy_sam_header, map_isoforms_and_sort, \
        can_merge, compare_fuzzy_junctions, collapse_fuzzy_junctions
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR

_SIV_DIR_ = op.join(SIV_DATA_DIR, "test_collapsing")
_DAT_DIR_ = op.join(DATA_DIR, "test_collapsing")
_OUT_DIR_ = op.join(OUT_DIR, "test_collapsing")

GMAP_INPUT_FASTA = op.join(_SIV_DIR_, 'gmap-input.fasta')
GMAP_INPUT_FASTQ = op.join(_SIV_DIR_, 'gmap-input.fastq')
GMAP_INPUT_FASTA_DS = op.join(_SIV_DIR_, 'gmap-input.fasta.contigset.xml')
GMAP_INPUT_FASTQ_DS = op.join(_SIV_DIR_, 'gmap-input.fastq.contigset.xml')
GMAP_SAM = op.join(SIV_DATA_DIR, 'test_SAMReader', 'gmap-output.sam')

GMAP_DB = op.join(SIV_DATA_DIR, 'gmap_db')
GMAP_NAME = 'SIRV'


class TEST_CollapsingUtils(unittest.TestCase):
    """Test functions of pbtranscript.collapsing.CollapsingUtils."""
    def setUp(self):
        """Define input and output file."""
        rmpath(_OUT_DIR_)
        mkdir(_OUT_DIR_)

    def test_copy_sam_header(self):
        """Test copy_sam_header"""
        out_fn = op.join(_OUT_DIR_, 'test_copy_sam_header.sam')
        copy_sam_header(GMAP_SAM, out_fn)
        with open(out_fn, 'r') as reader:
            lines = [r for r in reader]
            self.assertTrue(all([line.startswith('@') for line in lines]))
            self.assertEqual(len(lines), 9)

    def test_map_isoforms_and_sort(self):
        """Test map_isoforms_and_sort"""
        out_fn = op.join(_OUT_DIR_, 'test_map_isoforms_and_sort_fasta.sam')
        rmpath(out_fn)
        map_isoforms_and_sort(input_filename=GMAP_INPUT_FASTA,
                              sam_filename=out_fn,
                              gmap_db_dir=GMAP_DB,
                              gmap_db_name=GMAP_NAME,
                              gmap_nproc=10)
        self.assertTrue(op.exists(out_fn))

        out_fn = op.join(_OUT_DIR_, 'test_map_isoforms_and_sort_fastq.sam')
        rmpath(out_fn)
        map_isoforms_and_sort(input_filename=GMAP_INPUT_FASTQ,
                              sam_filename=out_fn,
                              gmap_db_dir=GMAP_DB,
                              gmap_db_name=GMAP_NAME,
                              gmap_nproc=10)
        self.assertTrue(op.exists(out_fn))

        out_fn = op.join(_OUT_DIR_, 'test_map_isoforms_and_sort_fasta_ds.sam')
        rmpath(out_fn)
        map_isoforms_and_sort(input_filename=GMAP_INPUT_FASTA_DS,
                              sam_filename=out_fn,
                              gmap_db_dir=GMAP_DB,
                              gmap_db_name=GMAP_NAME,
                              gmap_nproc=10)
        self.assertTrue(op.exists(out_fn))

        out_fn = op.join(_OUT_DIR_, 'test_map_isoforms_and_sort_fastq_ds.sam')
        rmpath(out_fn)
        map_isoforms_and_sort(input_filename=GMAP_INPUT_FASTQ_DS,
                              sam_filename=out_fn,
                              gmap_db_dir=GMAP_DB,
                              gmap_db_name=GMAP_NAME,
                              gmap_nproc=10)
        self.assertTrue(op.exists(out_fn))

    def test_collapse_fuzzy_junctions(self):
        """Test collapse_fuzzy_junctions, can_merge and compare_fuzzy_junctions."""
        test_name = "collapse_fuzzy_junctions"
        input_gff = op.join(_DAT_DIR_, "input_%s.gff" % test_name)
        input_group = op.join(_DAT_DIR_, "input_%s.group.txt" % test_name)
        output_gff = op.join(_OUT_DIR_, "output_%s.gff" % test_name)
        output_group = op.join(_OUT_DIR_, "output_%s.group.txt" % test_name)

        records = [r for r in CollapseGffReader(input_gff)]
        self.assertEqual(len(records), 4)

        r0, r1, r2, r3 = records
        # comparing r0 and r1
        m = compare_fuzzy_junctions(r0.ref_exons, r1.ref_exons, max_fuzzy_junction=5)
        self.assertEqual(m, "subset")
        self.assertTrue(can_merge(m, r0, r1, allow_extra_5exon=True, max_fuzzy_junction=5))

        # comparing r2 and r3
        m = compare_fuzzy_junctions(r2.ref_exons, r3.ref_exons, max_fuzzy_junction=5)
        self.assertEqual(m, "exact")
        self.assertTrue(can_merge(m, r2, r3, allow_extra_5exon=True, max_fuzzy_junction=5))

        # call collapse_fuzzy_junctions and write fuzzy output.
        collapse_fuzzy_junctions(gff_filename=input_gff,
                                 group_filename=input_group,
                                 fuzzy_gff_filename=output_gff,
                                 fuzzy_group_filename=output_group,
                                 allow_extra_5exon=True,
                                 max_fuzzy_junction=5)

        r4, r5 = [r for r in CollapseGffReader(output_gff)]
        self.assertEqual(r1, r4)
        self.assertEqual(r3, r5)
