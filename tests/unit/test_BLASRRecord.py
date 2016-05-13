"""Test classes defined within pbtranscript.io.BLASRRecord."""
import unittest
import os.path as op
import hashlib
from pbtranscript.io.BLASRRecord import BLASRM4Reader, \
    BLASRM5Reader, BLASRRecord
from test_setpath import DATA_DIR

class TEST_BLASRRECORD(unittest.TestCase):
    """Test classes defined within pbtranscript.io.BLASRRecord."""
    def setUp(self):
        """Define BLASR M4, M5 files."""
        self.dataDir = DATA_DIR
        self.m4 = op.join(self.dataDir, "test_BLASRRecord.m4")
        self.m5 = op.join(self.dataDir, "test_BLASRRecord.m5")
        self.t0 = BLASRRecord(
            qID=(
            "m130812_185809_42141_c100533960310000001823079711101380_s1_p0" +
            "/223/0_5667"),
            qLength=15693, qStart=10, qEnd=378, qStrand='+',
            sID="gi|49175990|ref|NC_000913.2|",
            sLength=4639675, sStart=3004911, sEnd=3005261, sStrand='-',
            score=-1489, mapQV=0, qAln=None, alnStr=None, sAln=None,
            strand='-', identity=89.8667)

        self.t1 = BLASRRecord(
            qID=(
            "m130812_185809_42141_c100533960310000001823079711101380_s1_p0" +
            "/223/0_5667"),
            qLength=15693, qStart=10, qEnd=392, qStrand='+',
            sID="gi|49175990|ref|NC_000913.2|",
            sLength=4639675, sStart=579683, sEnd=580047, sStrand='+',
            score=-1434, mapQV=254, qAln=None, alnStr=None, sAln=None,
            strand='+', identity=87.1795)

        self.t00 = BLASRRecord(
            qID=(
            "m130812_185809_42141_c100533960310000001823079711101380_s1_p0" +
            "/223/0_5667"),
            qLength=15693, qStart=10, qEnd=378, qStrand='+',
            sID="gi|49175990|ref|NC_000913.2|",
            sLength=4639675, sStart=4639675-3005261,
            sEnd=4639675-3004911, sStrand='-',
            score=-1489, mapQV=0,
            nMatch=337, nMismatch=6, nIns=25, nDel=7,
            strand='-', identity=89.8667)

        self.t11 = BLASRRecord(
            qID=(
            "m130812_185809_42141_c100533960310000001823079711101380_s1_p0" +
            "/223/0_5667"),
            qLength=15693, qStart=10, qEnd=392, qStrand='+',
            sID="gi|49175990|ref|NC_000913.2|",
            sLength=4639675, sStart=579683, sEnd=580047, sStrand='+',
            score=-1434, mapQV=254,
            nMatch=340, nMismatch=16, nIns=26, nDel=8,
            strand='+', identity=87.1795)

        self.t22 = BLASRRecord(
            qID=(
            "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/" +
            "293/3629_59_CCS"),
            qLength=3570, qStart=5, qEnd=3561, qStrand='+',
            sID=(
            "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/" +
            "27129/3686_54_CCS"),
            sLength=3632, sStart=5, sEnd=3621, sStrand='+',
            score=-15590, nMatch=3418, nMismatch=45, nIns=93, nDel=153,
            mapQV=0, qAln=None, alnStr=None, sAln=None,
            strand='+', identity=92.1542)

    def test_M4Reader(self):
        """Test BLASR M4 Reader."""
        print self.m4
        reader = BLASRM4Reader(self.m4)
        reads = [x for x in reader]
        self.assertTrue(len(reads), 2)
        r0, r1 = reads

        print "r0"
        print r0.__dict__
        print "t0"
        print self.t0.__dict__

        self.assertTrue(r0 == self.t0)
        self.assertTrue(r1 == self.t1)

    def test_M5Reader(self):
        """Test BLASR M5 Reader."""
        reader = BLASRM5Reader(self.m5)
        reads = [x for x in reader]
        self.assertTrue(len(reads), 3)

        r0, r1, r2 = reads
        self.assertTrue(hashlib.md5(r0.qAln).hexdigest() ==
            "7dfa7031ee5442f3b637248c0cb146dc")
        self.assertTrue(hashlib.md5(r0.alnStr).hexdigest() ==
            "8efa5e0c00482b66d4caaae0a530c791")
        self.assertTrue(hashlib.md5(r0.sAln).hexdigest() ==
            "4853820e603e032aa351ba14ae339032")

        r0.qAln = r0.alnStr = r0.sAln = None
        self.assertTrue(r0 == self.t00)

        self.assertTrue(hashlib.md5(r1.qAln).hexdigest() ==
            "efc9648fd043eb82be9c20f6a5f2a6dd")
        self.assertTrue(hashlib.md5(r1.alnStr).hexdigest() ==
            "7fa3f901b6173b998825fd75a8b201db")
        self.assertTrue(hashlib.md5(r1.sAln).hexdigest() ==
            "b12828acfb305b4e854b7d3b4da42010")

        r1.qAln = r1.alnStr = r1.sAln = None
#        r1.nMatch = r1.nMismatch = r1.nIns = r1.nDel = None
        self.assertTrue(r1 == self.t11)

        self.assertTrue(hashlib.md5(r2.qAln).hexdigest() ==
            "40d21d6f8160aca30c50c341ad37375a")
        self.assertTrue(hashlib.md5(r2.alnStr).hexdigest() ==
            "a7377bc9147d8b8297dc098df93145cc")
        self.assertTrue(hashlib.md5(r2.sAln).hexdigest() ==
            "8ce8f6a8055e8c54a76abd1eba893ea2")
        r2.qAln = r2.alnStr = r2.sAln = None
        print "r2"
        print r2.__dict__
        print "t22"
        print self.t22.__dict__
        self.assertTrue(r2 == self.t22)


