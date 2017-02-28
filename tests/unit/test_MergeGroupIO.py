#!/usr/bin/env python
"""
Test MegaInfoReader/Writer
"""

import unittest
import os
import os.path as op

from pbtranscript.io import MegaInfoReader, MegaInfoWriter, MergeGroupOperation
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR


MEGA_FN = op.join(DATA_DIR, 'test.mega_info.txt')


class TestMegaInfoIO(unittest.TestCase):

    def test_all(self):
        """Test MegaInfoReader, Writer"""
        cwd = os.getcwd()
        rs = [r for r in MegaInfoReader(MEGA_FN)]
        self.assertEqual([r['pbid'] for r in rs], ['PB.1.1', 'PB.1.2', 'PB.2.1', 'PB.3.1'])
        self.assertEqual([r['SF3BI'] for r in rs], ['NA', 'PB.1.1', 'NA', 'PB.4.1'])
        self.assertEqual([r['NTI'] for r in rs], ['PB.2.1', 'NA', 'PB.5.1', 'PB.7.1'])

        o_mega_fn = op.join(OUT_DIR, 'o.mega_info.txt')
        writer = MegaInfoWriter(o_mega_fn, 'SF3BI', 'NTI')
        writer.writeRecord(MergeGroupOperation('PB.1.1', None, 'PB.2.1'))
        writer.writeRecord(MergeGroupOperation('PB.1.2', 'PB.1.1', None))
        writer.writeRecord(MergeGroupOperation('PB.2.1', None, 'PB.5.1'))
        writer.writeRecord(MergeGroupOperation('PB.3.1', 'PB.4.1', 'PB.7.1'))

        writer.close()
        self.assertTrue(op.exists(o_mega_fn))

        ors = [r for r in MegaInfoReader(o_mega_fn)]
        self.assertEqual(str(rs), str(ors))
