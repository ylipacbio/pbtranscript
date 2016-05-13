#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbtranscript.ClusterOptions import SgeOptions

class TestSgeOptions(unittest.TestCase):
    """Test pbtranscript.ClusterOptions.SgeOptions"""
    def setUp(self):
        """Initialize."""
        pass

    def test_qsub_cmd(self):
        """Test qsub_cmd."""
        sge_opts = SgeOptions(unique_id=100)
        self.assertEqual(sge_opts.qsub_cmd("a.sh", num_threads=1),
                         "qsub -cwd -V -S /bin/bash -pe smp 1 -e /dev/null -o /dev/null a.sh")

        sge_opts = SgeOptions(unique_id=100, sge_queue="my_sge_queue",
                              sge_env_name="orte")

        self.assertEqual(sge_opts.qsub_cmd("a.sh", num_threads=1,
                         wait_before_exit=True, depend_on_jobs=['1', '2', '3']),
                         "qsub -cwd -V -S /bin/bash -pe orte 1 -q my_sge_queue -sync y -hold_jid 1,2,3 -e /dev/null -o /dev/null a.sh")
