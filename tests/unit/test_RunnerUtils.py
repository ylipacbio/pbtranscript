#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbcore.util.Process import backticks
from pbtranscript.RunnerUtils import *
from pbtranscript.ClusterOptions import SgeOptions
from test_setpath import DATA_DIR, OUT_DIR, STD_DIR, SIV_DATA_DIR

class TestRunnerUtils(unittest.TestCase):
    """Test pbtranscript.RunnerUtils"""
    def setUp(self):
        """Initialize."""
        self.data_dir = DATA_DIR
        self.out_dir = OUT_DIR
        self.stdout_dir = STD_DIR
        self.sivDataDir = SIV_DATA_DIR

    def test_write_cmd_to_script(self):
        """Test write_cmd_to_script."""
        expected_content = ["#!/bin/bash\n", "echo hello\n"]

        outfn = op.join(self.out_dir, "test_write_cmd_to_script.sh")
        write_cmd_to_script("echo hello", outfn)
        content = [r for r in open(outfn, 'r')]
        self.assertEqual(content, expected_content)

    def test_local_job_runner(self):
        """Test local_job_runner."""
        cmds_list = ["echo 1", "echo 2", "echo 3"]
        num_threads = 2
        local_job_runner(cmds_list, num_threads)

        cmds_list.append("unknown_cmd")
        detectError = False

        try:
            local_job_runner(cmds_list, num_threads, throw_error=True)
        except RuntimeError:
            detectError = True

        self.assertTrue(detectError)

        self.assertEqual(["unknown_cmd"],
                         local_job_runner(cmds_list, num_threads, throw_error=False))

    @unittest.skipUnless(backticks('qstat')[1] == 0, "sge disabled")
    def test_get_active_sge_jobs(self):
        """Test get_active_sge_jobs"""
        self.assertTrue(isinstance(get_active_sge_jobs(), dict))

    @unittest.skipUnless(backticks('qstat')[1] == 0, "sge disabled")
    def test_sge_submit(self):
        """Test sge_submit."""
        jid = sge_submit("qsub -cwd -b y -e /dev/null -o /dev/null echo 1")
        self.assertTrue(isinstance(jid, str))
        self.assertTrue(isinstance(int(jid), int))
        # no need to wait for it to finish

    @unittest.skipUnless(backticks('qstat')[1] == 0, "sge disabled")
    def test_sge_job_runner(self):
        """Test sge_job_runner"""
        cmds = ["sleep 5", "sleep 5", "sleep 5", "sleep 5"]
        script_files = [op.join(self.out_dir,
                                "test_sge_job_runner_%s.sh" % i)
                        for i in range(0, len(cmds))]
        #done_script = op.join(self.out_dir, "test_sge_job_runner.done.sh")
        #done_file = op.join(self.out_dir, "test_sge_job_runner.done")

        delete_files = script_files #+ [done_script, done_file]
        for f in delete_files:
            backticks('rm %s' % f)

        sge_opts = SgeOptions(100)
        #write_cmd_to_script(cmd="echo 'done' > %s" % done_file,
        #                    script=done_script)
        jids = sge_job_runner(cmds, script_files=script_files,
                              #done_script=done_script,
                              num_threads_per_job=1,
                              sge_opts=sge_opts, qsub_try_times=1)
        #self.assertTrue(op.exists(done_file))

