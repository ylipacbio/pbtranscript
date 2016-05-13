#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbtranscript.ClusterOptions import IceOptions
from pbtranscript.Utils import mknewdir, execute
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR, SIV_STD_DIR


def copy_in_fasta_to_out(in_dir, out_dir, filename):
    """copy filename from in_dir (e.g., data) to out_dir,
    return out_fasta
    """
    mknewdir(out_dir)
    cmd = "cp %s %s" % (op.join(in_dir, filename),
                        op.join(out_dir, filename))
    execute(cmd=cmd)
    return op.join(out_dir, filename)


class TestIceOptions(unittest.TestCase):
    """Test pbtranscript.ClusterOptions.IceOptions"""
    def setUp(self):
        """Initialize."""
        self.filename = "reads_of_insert.fasta"

    def test_read_write_config(self):
        """test read_config and write_config."""
        out_dir = op.join(OUT_DIR, "test_ice_opts_read_write_config")
        fasta_filename = copy_in_fasta_to_out(DATA_DIR, out_dir, self.filename)

        ice_opts = IceOptions(cDNA_size="above5k")
        ice_opts._write_config(fasta_filename)
        self.assertTrue(op.exists(ice_opts._config_filename(
            fasta_filename=fasta_filename)))

        self.assertEqual(ice_opts.low_cDNA_size, 739)
        self.assertEqual(ice_opts.high_cDNA_size, 4175)
        self.assertEqual(ice_opts.sensitive_mode, False)

        ice_opts._read_config(fasta_filename)

        self.assertEqual(ice_opts.low_cDNA_size, 739)
        self.assertEqual(ice_opts.high_cDNA_size, 4175)
        self.assertEqual(ice_opts.sensitive_mode, False)

    def test_detect_cDNA_size(self):
        """Test detect_cDNA_size."""
        out_dir = op.join(OUT_DIR, "test_ice_opts_detect_cDNA_size")
        fasta_filename = copy_in_fasta_to_out(DATA_DIR, out_dir, self.filename)

        ice_opts = IceOptions(cDNA_size="above5k")

        config_filename = ice_opts._config_filename(fasta_filename=fasta_filename)
        if op.exists(config_filename):
            os.remove(config_filename)

        ice_opts.detect_cDNA_size(fasta_filename)

        self.assertTrue(op.exists(config_filename))

        self.assertEqual(ice_opts.low_cDNA_size, 739)
        self.assertEqual(ice_opts.high_cDNA_size, 4175)
        self.assertEqual(ice_opts.cDNA_size, "under1k")
        self.assertEqual(ice_opts.sensitive_mode, False)
