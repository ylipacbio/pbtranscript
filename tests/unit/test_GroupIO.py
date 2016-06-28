"""Test pbtranscript.branch."""
import unittest
import os.path as op
import filecmp
from pbtranscript.io.GroupIO import GroupRecord, GroupReader, GroupWriter
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR, SIV_STD_DIR

GROUP_FN_1 = op.join(DATA_DIR, 'test_GroupReader.txt')
GROUP_FN_2 = op.join(SIV_STD_DIR, 'test_branch/test_branch.group.txt.fuzzy')

class TEST_GroupIO(unittest.TestCase):
    """Test pbtranscript.io.GroupRecord, GroupReader, GroupWriter."""
    def setUp(self):
        """Define input and output file."""
        pass

    def test_all(self):
        """Test All"""
        expected_r = GroupRecord(name="group1",
                                 members=["member0", "member1", "member2"])
        with GroupReader(GROUP_FN_1) as reader:
            records = [r for r in reader]
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0], expected_r)

        expected_r = GroupRecord(name="PB.1.1",
                                 members="i0_HQ_sampleb92221|c8319/f2p0/463,i0_HQ_sampleb92221|c28/f4p0/460,i0_HQ_sampleb92221|c524/f2p0/462,i0_HQ_sampleb92221|c539/f2p0/460,i0_HQ_sampleb92221|c7864/f22p0/462,i0_HQ_sampleb92221|c7959/f2p0/461,i0_HQ_sampleb92221|c8090/f3p0/462,i0_HQ_sampleb92221|c8099/f3p0/459,i0_HQ_sampleb92221|c8136/f2p0/461,i0_HQ_sampleb92221|c428/f2p0/459".split(','))
        with GroupReader(GROUP_FN_2) as reader:
            records = [r for r in reader]
            self.assertEqual(len(records), 51)
            self.assertEqual(records[0], expected_r)
            out_fn = op.join(OUT_DIR, "test_GroupWriter.txt")
            with GroupWriter(out_fn) as writer:
                for r in records:
                    writer.writeRecord(r)

            self.assertTrue(filecmp.cmp(out_fn, GROUP_FN_2))
