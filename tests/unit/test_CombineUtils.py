"""Test classes defined within pbtranscript.CombineUtils."""

import unittest
import os.path as op

from pbcore.io import FastaReader, FastqReader
from pbtranscript.Utils import rmpath, mkdir
from pbtranscript.ClusterOptions import IceQuiverHQLQOptions
from pbtranscript.CombineUtils import CombineRunner
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR


class TEST_CombineUtils(unittest.TestCase):

    """Test functions of pbtranscript.CombineUtils."""

    def setUp(self):
        """Define input and output file."""
        pass

    def test_runner(self):
        """Test CombineRunner."""
        ipq_opts = IceQuiverHQLQOptions(qv_trim_5=100, qv_trim_3=30)
        d = op.join(SIV_DATA_DIR, "test_tool_contract_chunks")
        split_dirs = [op.join(d, b, "cluster_out") for b in
                      ("0to1kb_part0", "1to2kb_part0", "2to3kb_part0", "3to4kb_part0", "4to5kb_part0")]
        print split_dirs
        out_combined_dir = op.join(OUT_DIR, "test_CombineUtils", "combined_dir")
        rmpath(out_combined_dir)
        mkdir(out_combined_dir)
        obj = CombineRunner(combined_dir=out_combined_dir,
                            sample_name="mysample",
                            split_dirs=split_dirs,
                            ipq_opts=ipq_opts)
        obj.run()

        expected_out_fns = (obj.all_hq_fa, obj.all_hq_fq, obj.all_lq_fa, obj.all_lq_fq,
                            obj.all_consensus_isoforms_fa,
                            obj.all_cluster_report_fn, obj.all_cluster_summary_fn)
        self.assertTrue(all([op.exists(f) for f in expected_out_fns]))

        expected_hq_isoforms = ['i1_HQ_mysample|c0/f2p16/1826', 'i2_HQ_mysample|c2/f9p14/2470',
                                'i2_HQ_mysample|c5/f7p19/2472', 'i2_HQ_mysample|c10/f8p16/2457',
                                'i2_HQ_mysample|c98/f2p10/2081', 'i2_HQ_mysample|c108/f23p28/2471']
        self.assertEqual([r.name.split(' ')[0] for r in FastaReader(obj.all_hq_fa)], expected_hq_isoforms)
        self.assertEqual([r.name.split(' ')[0] for r in FastqReader(obj.all_hq_fq)], expected_hq_isoforms)

        expected_lq_isoforms_num = 73
        self.assertEqual(len([r for r in FastaReader(obj.all_lq_fa)]), expected_lq_isoforms_num)

        expected_consensus_isoforms_num = 79
        self.assertEqual(len([r for r in FastaReader(obj.all_consensus_isoforms_fa)]), expected_consensus_isoforms_num)

