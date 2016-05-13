#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
import shutil
from pbtranscript.Utils import cat_files, filter_sam, validate_fofn, \
        get_sample_name, mknewdir, as_contigset, execute
from test_setpath import DATA_DIR, OUT_DIR, STD_DIR, SIV_DATA_DIR

class TestUtils(unittest.TestCase):
    """Test pbtranscript.Utils"""
    def setUp(self):
        """Initialize."""
        self.data_dir = DATA_DIR
        self.out_dir = OUT_DIR
        self.stdout_dir = STD_DIR
        self.sivDataDir = SIV_DATA_DIR

    def test_cat_files(self):
        """Test cat_files."""
        fn_1 = op.join(self.data_dir, "primers.fasta")
        fn_2 = op.join(self.data_dir, "test_phmmer.fasta")
        out_fn_1 = op.join(self.out_dir, "test_cat_1")
        out_fn_2 = op.join(self.out_dir, "test_cat_2")

        std_out_fn_2 = op.join(self.stdout_dir, "test_cat_2")

        cat_files(src=[fn_1], dst=out_fn_1)
        cat_files(src=[fn_1, fn_2], dst=out_fn_2)
        self.assertTrue(filecmp.cmp(out_fn_1, fn_1))
        self.assertTrue(filecmp.cmp(out_fn_2, std_out_fn_2))

    def test_filter_sam(self):
        """Test filter_sam."""
        in_sam = op.join(self.data_dir, "test_filter_sam.sam")
        out_sam = op.join(self.out_dir, "test_filter_sam.sam")
        stdout_sam = op.join(self.stdout_dir, "test_filter_sam.sam")
        filter_sam(in_sam=in_sam, out_sam=out_sam)
        self.assertTrue(filecmp.cmp(out_sam, stdout_sam))

    def test_validate_fofn(self):
        """Test validate_fofn"""
        detectedIOError = False
        fofn_filename = op.join(self.data_dir, "bad.fofn")
        try:
            validate_fofn(fofn_filename)
        except IOError as e:
            detectedIOError = True
        self.assertTrue(detectedIOError)

    def test_guess_file_format(self):
        """Test guess_file_format."""
        from pbtranscript.Utils import guess_file_format, FILE_FORMATS
        fn1 = op.join(self.sivDataDir, "bigbam", "ccsbam.fofn")
        self.assertTrue(guess_file_format(fn1) == FILE_FORMATS.BAM)

        fn2 = op.join(self.sivDataDir, "bigbam", "bas.fofn")
        self.assertTrue(guess_file_format(fn2) == FILE_FORMATS.H5)

        self.assertTrue(guess_file_format([fn1, fn2]) == FILE_FORMATS.UNKNOWN)

    def test_get_sample_name(self):
        """Test get_sample_name"""
        self.assertTrue(get_sample_name("my_name"), "my_name")
        self.assertTrue(get_sample_name("my name,|"), "myname")
        self.assertTrue(len(get_sample_name("")) > 0)

    def test_as_contigset(self):
        """Test as_contigset"""
        out_dir = op.join(OUT_DIR, 'test_Utils')
        mknewdir(out_dir)
        fa = op.join(out_dir, "empty.fasta")
        xml = op.join(out_dir, "empty.contigset.xml")
        fai = fa + ".fai"

        execute("touch %s" % fa)
        as_contigset(fa, xml)
        self.assertTrue(op.exists(xml))
        self.assertTrue(op.exists(fai))

        fn = 'reads_of_insert.fasta'
        shutil.copy(src=op.join(DATA_DIR, fn), dst=op.join(out_dir, fn))
        fa = op.join(out_dir, fn)
        as_contigset(fa, fa)

        fai = fa + ".fai"
        xml = op.join(out_dir, 'reads_of_insert.contigset.xml')
        as_contigset(fa, xml)
        self.assertTrue(op.exists(xml))
        self.assertTrue(op.exists(fai))

