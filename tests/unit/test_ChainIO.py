#!/usr/bin/env python
"""
Test ChainConfig
"""

import unittest
import os
import os.path as op

from pbtranscript.io.ChainIO import ChainConfig, SampleFiles
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR


CFG_FN = op.join(DATA_DIR, 'sample.config')


class TestChainConfig(unittest.TestCase):

    def test_all(self):
        """Test ChainConfig"""
        cwd = os.getcwd()
        group_fn =  'hq_isoforms.collapsed.group.txt'
        gff_fn = 'hq_isoforms.collapsed.min_fl_2.filtered.gff'
        abundance_fn = 'hq_isoforms.collapsed.min_fl_2.filtered.abundance.txt'
        samples_dict = {'SF3BI': op.join(cwd, 'SF3BI'), 'NTI': op.join(cwd, 'NTI')}

        cfg = ChainConfig.from_file(CFG_FN)
        self.assertTrue(len(cfg.sample_names) == 2)
        for name, path in zip(cfg.sample_names, cfg.sample_paths):
            self.assertTrue(name in samples_dict)
            self.assertEqual(path, samples_dict[name])

        self.assertEqual(cfg.group_fn, group_fn)
        self.assertEqual(cfg.gff_fn, gff_fn)
        self.assertEqual(cfg.abundance_fn, abundance_fn)

        path = op.join(cwd, 'NTI')
        sample0 = SampleFiles(name='NTI', path=path, group_fn=op.join(path, group_fn), gff_fn=op.join(path, gff_fn), abundance_fn=op.join(path, abundance_fn))
        self.assertEqual(str(cfg.samples[0].name), 'SF3BI') # the first sample must be SF3BI
        self.assertEqual(str(cfg.samples[1]), str(sample0))

        o_cfg_fn = op.join(OUT_DIR, 'sample.config')
        cfg.write(o_cfg_fn)
        self.assertTrue(op.exists(o_cfg_fn))
        ocfg = ChainConfig.from_file(o_cfg_fn)
        self.assertEqual(str(ocfg), str(cfg))
