#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbcore.util.Process import backticks
from pbtranscript.ClusterOptions import SgeOptions
from pbtranscript.Utils import mkdir, mknewdir
from pbtranscript.ice_daligner import DalignerRunner
from test_setpath import DATA_DIR, OUT_DIR, STD_DIR, SIV_DATA_DIR

class TestDalignerRunner(unittest.TestCase):
    """Test pbtranscript.ice_daligner.DalignerRunner"""
    def setUp(self):
        """Initialize."""
        self.data_dir = op.join(DATA_DIR, "test_daligner_against_ref")
        self.script_dir = op.join(OUT_DIR, "test_ice_daligner_script")
        self.dazz_dir = op.join(OUT_DIR, "test_ice_daligner_dazz")
        self.out_dir = op.join(OUT_DIR, "test_ice_daligner_out")
        mkdir (self.dazz_dir)
        mkdir (self.out_dir)
        self.stdout_dir = STD_DIR
        self.sivDataDir = SIV_DATA_DIR
        self.query_filename  = "test_daligner_query.fasta"
        self.target_filename = "test_daligner_target.fasta"
        self.runner = DalignerRunner(query_filename=op.join(self.data_dir, self.query_filename),
                                     target_filename=op.join(self.data_dir, self.target_filename),
                                     is_FL=False, same_strand_only=True,
                                     dazz_dir=self.dazz_dir, script_dir=self.script_dir)
        self.runner.output_dir = self.out_dir

    def test_query_prefix(self):
        """Test query_prefix."""
        self.assertEqual(self.runner.query_prefix(1),
                op.join(self.dazz_dir, self.query_filename[0:-6] + ".dazz.fasta"))

    def test_target_prefix(self):
        """Test target_prefix."""
        self.assertEqual(self.runner.target_prefix(1),
                op.join(self.dazz_dir, self.target_filename[0:-6] + ".dazz.fasta"))

    def test_thread_prefix(self):
        """Test local_job_runner."""
        self.assertEqual(self.runner.thread_prefix(2, is_forward=True), 'N2')
        self.assertEqual(self.runner.thread_prefix(2, is_forward=False), 'C2')

    def test_las_filenames(self):
        """Test las_filenames."""
        expected = [op.join(self.out_dir,
                            "{q}.{t}.{k}.las".format(
                                q=self.query_filename[0:-6] + ".dazz.fasta",
                                t=self.target_filename[0:-6] + ".dazz.fasta",
                                k=k))
                    for k in ('N0', 'N1', 'N2', 'N3')]
        print 'las_filenames\n'
        print expected
        self.assertEqual(self.runner.las_filenames, expected)

    def test_la4ice_filenames(self):
        """Test la4ice_filenames."""
        expected = [op.join(self.out_dir,
                            "{q}.{t}.{k}.las.out".format(
                                q=self.query_filename[0:-6] + ".dazz.fasta",
                                t=self.target_filename[0:-6] + ".dazz.fasta",
                                k=k))
                    for k in ('N0', 'N1', 'N2', 'N3')]
        self.assertEqual(self.runner.la4ice_filenames, expected)

    def test_run(self):
        """Test run(output_dir, min_match_len, sensitive_mode).
        running on sge and locally.
        """
        run_on_sge = (backticks('qstat')[1] == 0)

        if run_on_sge:
            self.runner.use_sge = True
            self.runner.sge_opts = SgeOptions(100)
            mknewdir(self.out_dir)
            self.runner.run(output_dir=self.out_dir)

            for las_filename in self.runner.las_filenames:
                print "Checking existance of " + las_filename
                self.assertTrue(op.exists(las_filename))

            for la4ice_filename in self.runner.la4ice_filenames:
                print "Checking existance of " + la4ice_filename
                self.assertTrue(op.exists(la4ice_filename))

        # Run locally
        self.runner.use_sge = False
        mknewdir(self.out_dir)
        self.runner.run(output_dir=self.out_dir)

        for las_filename in self.runner.las_filenames:
            print "Checking existance of " + las_filename
            self.assertTrue(op.exists(las_filename))

        for la4ice_filename in self.runner.la4ice_filenames:
            print "Checking existance of " + la4ice_filename
            self.assertTrue(op.exists(la4ice_filename))

        # clean all output
        self.runner.clean_run()

        for las_filename in self.runner.las_filenames:
            print "Checking %s has been removed.\n" % las_filename
            self.assertTrue(not op.exists(las_filename))

        for la4ice_filename in self.runner.la4ice_filenames:
            print "Checking %s has been removed.\n" % la4ice_filename
            self.assertTrue(not op.exists(la4ice_filename))
