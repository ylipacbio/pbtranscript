"""Test pbtranscript.Classifier."""

import unittest
import filecmp
import os.path as op

from pbcommand.pb_io.report import load_report_from_json
from pbcore.util.Process import backticks

from pbtranscript.io.Summary import ClassifySummary, ClusterSummary


def rm_version_string(infn, outfn):
    cmd = "cat %s |grep -v '_changelist' | grep -v '_version' > %s" % (infn, outfn)
    _o, _c, _e = backticks(cmd)
    if _c != 0:
        raise RuntimeError("Failed to run %s" % cmd)

def _compare_reports(self, rpt_json1, rpt_json2):
    rpt1 = load_report_from_json(rpt_json1)
    rpt2 = load_report_from_json(rpt_json2)
    attr1 = {a.id:a.value for a in rpt1.attributes}
    attr2 = {a.id:a.value for a in rpt2.attributes}
    self.assertEqual(attr1, attr2)

class Test_ClassifySummary(unittest.TestCase):
    """Test ClassifySummary."""
    def setUp(self):
        """Set up test data."""
        self.rootDir = op.dirname(op.dirname(op.abspath(__file__)))
        self.testDir = op.join(self.rootDir, "")

    def test_write(self):
        """Test ClassifySummary.write."""
        obj = ClassifySummary()
        obj.num_reads = 100
        obj.num_5_seen = 90
        obj.num_3_seen = 70
        obj.num_polya_seen = 70
        obj.num_filtered_short_reads = 10
        obj.num_nfl = 50
        obj.num_fl = 40
        obj.num_flnc = 39
        obj.num_flc = 1
        obj.num_flnc_bases = 10001

        outFN = op.join(self.testDir, "out/test_ClassifySummary.txt")
        stdoutFN = op.join(self.testDir, "stdout/test_ClassifySummary.txt")
        obj.write(outFN)
        self.assertTrue(filecmp.cmp(outFN, stdoutFN))

        outFN = op.join(self.testDir, "out/test_ClassifySummary.json")
        stdoutFN = op.join(self.testDir, "stdout/test_ClassifySummary.json")
        obj.write(outFN)
        rm_version_string(outFN, outFN + "tmp1")
        rm_version_string(stdoutFN, outFN + "tmp2")
        _compare_reports(self, outFN, stdoutFN)
        #self.assertTrue(filecmp.cmp(outFN + "tmp1", outFN + "tmp2"))


class Test_ClusterSummary(unittest.TestCase):
    """Test ClusterSummary."""
    def setUp(self):
        """Set up test data"""
        self.testDir = op.dirname(op.dirname(op.abspath(__file__)))

    def test_write(self):
        """Test ClusterSummary.write."""
        obj = ClusterSummary()
        obj.num_consensus_isoforms = 97
        obj.num_total_bases = 97 * 3945

        outFN = op.join(self.testDir, "out/test_ClusterSummary.txt")
        stdoutFN = op.join(self.testDir, "stdout/test_ClusterSummary.txt")
        obj.write(outFN)
        self.assertTrue(filecmp.cmp(outFN, stdoutFN))

        outFN = op.join(self.testDir, "out/test_ClusterSummary.json")
        stdoutFN = op.join(self.testDir, "stdout/test_ClusterSummary.json")
        obj.write(outFN)

        rm_version_string(outFN, outFN + "tmp1")
        rm_version_string(stdoutFN, outFN + "tmp2")
        _compare_reports(self, outFN, stdoutFN)
        #self.assertTrue(filecmp.cmp(outFN + "tmp1", outFN + "tmp2"))
