"""Test classes defined within pbtranscript.counting.CountUtils."""
import unittest
import os.path as op
from pbtranscript.Utils import rmpath, mkdir
from pbtranscript.counting.CountUtils import get_roi_len, read_group_filename
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR

_SIV_DIR_ = op.join(SIV_DATA_DIR, "test_counting")
_DAT_DIR_ = op.join(DATA_DIR, "test_counting")
_OUT_DIR_ = op.join(OUT_DIR, "test_counting")

GROUP_FN = op.join(_SIV_DIR_, "group.txt")

class TEST_CountUtils(unittest.TestCase):
    """Test functions of pbtranscript.counting.CountUtils."""
    def setUp(self):
        """Define input and output file."""
        rmpath(_OUT_DIR_)
        mkdir(_OUT_DIR_)

    def test_get_roi_len(self):
        """Test get_roi_len"""
        self.assertTrue(get_roi_len("movie/124/100_300_CCS"), 200)
        self.assertTrue(get_roi_len("movie/124/500_200_CCS"), 300)

    def test_read_group_filename(self):
        """Test read group filename."""
        cid_info = read_group_filename(group_filename=GROUP_FN,
                                       is_cid=True,
                                       sample_prefixes=None)

        self.assertEqual(len(cid_info[None].keys()), 846)
        self.assertEqual(cid_info[None]['i1_HQ_sampleb92221|c1030'], 'PB.5.6')
        self.assertEqual(cid_info[None]['i2_HQ_sampleb92221|c326'], 'PB.10.14')


        sample_prefixes = ['i0_HQ_sampleb92221', 'i1_HQ_sampleb92221', 'i2_HQ_sampleb92221']
        cid_info = read_group_filename(group_filename=GROUP_FN,
                                       is_cid=True,
                                       sample_prefixes=sample_prefixes)
        self.assertEqual(len(cid_info[sample_prefixes[0]].keys()), 165)
        self.assertEqual(len(cid_info[sample_prefixes[1]].keys()), 297)
        self.assertEqual(len(cid_info[sample_prefixes[2]].keys()), 384)
        self.assertEqual(cid_info['i1_HQ_sampleb92221']['c1030'], 'PB.5.6')
        self.assertEqual(cid_info['i2_HQ_sampleb92221']['c326'], 'PB.10.14')
