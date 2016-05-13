"""Test classes defined within pbtranscript.io.PbiBamIO."""
import unittest
import os.path as op
import hashlib
import time

from pbtranscript.Utils import make_pbi
from pbtranscript.io.PbiBamIO import BamCollection, BamHeader, BamWriter, BamZmwRead
from pbcore.io import BamAlignment, IndexedBamReader, readFofn, ConsensusReadSet
from test_setpath import OUT_DIR, STD_DIR
from pbcore.util.Process import backticks

def compareBamRecords(this, other):
    """Compare this (a BamAlignment object) with other
      (a PbiBamReader.BamZmwRead object)"""
    assert(isinstance(this, BamAlignment) and
           isinstance(other, BamZmwRead))
    return (this.readName == other.readName    and
            this.zmwName  == other.zmw.zmwName and
            this.aStart   == other.readStart   and
            this.aEnd     == other.readEnd     and
            this.read(False) == other.basecalls() and
            all(this.DeletionQV(False) == other.qv("DeletionQV")) and
            all(this.InsertionQV(False) == other.qv("InsertionQV")) and
            all(this.SubstitutionQV(False) == other.qv("SubstitutionQV"))
           )


def _verify_write_compare_subreads(testobj, inbamfns, zmws, outbamfn):
    """First verify that input bam and pbi files exist,
    next extract zmws from inputs and write to outbamfn,
    then compare bam records in input and output."""
    # Verify that input.bam and input.bam.pbi exist
    testobj.assertTrue(all(op.exists(fn) for fn in inbamfns))
    testobj.assertTrue(all(op.exists(fn + ".pbi") for fn in inbamfns))

    reader = BamCollection(*inbamfns)
    writer = BamWriter(outbamfn, reader.header)
    for zmw in zmws:
        for sr in reader[zmw].subreads:
            writer.write(sr)
    writer.close()

    # make pbi for outbamfn
    make_pbi(outbamfn)
    testobj.assertTrue(op.exists(outbamfn + ".pbi"))

    # Read subreads from outbamfn and compare.
    reader2 = IndexedBamReader(outbamfn)
    for r in reader2:
        other = reader[r.readName]
        testobj.assertTrue(compareBamRecords(r, other))


def _verify_write_compare_ccs(testobj, inbamfns, zmws, outbamfn,
                              expected_movies, expected_len):
    """First verify input.bam and input.bam.pbi exist,
    next, extract zmws from input and write to an output bam,
    then compare ccs reads and zmws in input and output.
    """
    testobj.assertTrue(all(op.exists(fn) for fn in inbamfns))
    testobj.assertTrue(all(op.exists(fn + ".pbi") for fn in inbamfns))

    reader = BamCollection(*inbamfns)
    # verify movie names and length of reader.

    testobj.assertTrue(set(reader.movieNames) == set(expected_movies))
    testobj.assertTrue(len(reader) == expected_len, "%d != %d" %(len(reader),
            expected_len))

    # write ccs reads.
    with BamWriter(outbamfn, reader.header) as writer:
        for zmw in zmws:
            writer.write(reader[zmw].ccsRead)

    # make pbi and check
    make_pbi(outbamfn)
    testobj.assertTrue(op.exists(outbamfn + ".pbi"))

    # compare ccs reads in input and output.
    reader2 = IndexedBamReader(outbamfn)
    outzmws = []
    for r in reader2:
        outzmws.append(r.zmwName)
        other = reader[r.readName]
        testobj.assertTrue(compareBamRecords(r, other))

    # compare ccs zmws in input and output
    testobj.assertTrue(set(zmws) == set(outzmws))


class TestIO(unittest.TestCase):
    def setUp(self):
        """Define input bam file, out bam file."""
        pDataDir = "/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data"
        self.dataDir = op.join(pDataDir, "bam")
        self.bigDataDir = op.join(pDataDir, "bigbam")
        self.ioDataDir = op.join(pDataDir, "io")
        self.outDir  = OUT_DIR
        self.stdDir  = STD_DIR

    def test_fofn(self):
        fofn = op.join(self.dataDir, "bam.fofn")
        bc = BamCollection(fofn) # XXX just making sure this works

    def test_read_subreads_from_one_file_of_one_smrtcell_write_bam(self):
        """
        test_read_subreads_from_one_file_of_one_smrtcell_write_bam
        Read subreads from a single bam file of a SMRTCell,
        select records of a few zmws and write to an output bam file.
        """
        movie = "m131018_081703_42161_c100585152550000001823088404281404_s1_p0"

        subreadsfns  = [op.join(self.dataDir, "%s.1.subreads.bam" % movie)]
        outbamfn = op.join(self.outDir, "test_writebam_subreads_1.bam")

        # zmws to extract
        hns = [45, 161, 227, 293, 495, 642, 780, 865, 888]
        zmws = ["%s/%s" % (movie, hn) for hn in hns]

        _verify_write_compare_subreads(self, subreadsfns, zmws, outbamfn)

    def test_read_subreads_from_multiple_files_of_one_smrtcell_write_bam(self):
        """
        test_read_subreads_from_multiple_files_of_one_smrtcell_write_bam
        Read raw subreads from multiple bam files of ONE SMRTCell,
        select records of a few zmws from each bam file and
        write to an output bam file.
        """
        movie = "m131018_081703_42161_c100585152550000001823088404281404_s1_p0"

        subreadsfns = [op.join(self.dataDir, "%s.%d.subreads.bam" % (movie, n)) for n in [1,2,3]]
        outbamfn = op.join(self.outDir, "test_writebam_subreads_2.bam")

        # zmws to extract
        hns = [45, 161, 227, 293, 495, 642, 780, 865, 888,
               54501, 56734,
               109028, 109075, 111128 ]
        zmws = ["%s/%s" % (movie, hn) for hn in hns]

        _verify_write_compare_subreads(self, subreadsfns, zmws, outbamfn)

    def test_read_subreads_from_multiple_files_of_multiple_smrtcells_write_bam(self):
        """
        test_read_subreads_from_multiple_files_of_multiple_smrtcells_write_bam
        Read subreads from multiple bam files of multiple SMRTCells,
        select records of a few zmws from each file and
        write to an output bam file.
        """
        movies = ["m140802_021814_42161_c110036822550000001823106706241500_s1_p0",
                  "m140802_043938_42161_c110036822550000001823106706241501_s1_p0",
                  "m140802_070303_42161_c110036822550000001823106706241502_s1_p0"]

        subreadsfns = [op.join(self.bigDataDir, "%s.1.subreads.bam" % movie) for movie in movies]
        outbamfn = op.join(self.outDir, "test_writebam_subreads_3.bam")

        hns = ["27",
               "1359",
               "988"]
        zmws = ["%s/%s" % (movie, hn) for (movie, hn) in zip(movies, hns)]

        _verify_write_compare_subreads(self, subreadsfns, zmws, outbamfn)

    def test_read_ccs_from_one_file_of_one_smrtcell_write_bam(self):
        """
        test_read_ccs_from_one_file_of_one_smrtcell_write_bam:
        Read ccs reads from one bam file of one SMRTCell,
        select records of a few zmws from each file, and
        write to an output bam file.
        """
        movies = ["m140802_021814_42161_c110036822550000001823106706241500_s1_p0"]

        ccsbamfns = [op.join(self.bigDataDir, "%s.1.ccs.bam" % movie) for movie in movies]
        outbamfn = op.join(self.outDir, "test_writebam_ccs_1.bam")

        hns = [58, 8473]
        zmws = ["%s/%s" % (movies[0], hn) for hn in hns]
        expected_num_ccs_reads = 13821

        _verify_write_compare_ccs(self, ccsbamfns, zmws, outbamfn,
                                  expected_movies=movies,
                                  expected_len=expected_num_ccs_reads)

    def test_read_ccs_from_multiple_files_of_one_smrtcell_write_bam(self):
        """
        test_read_ccs_from_multiple_files_of_one_smrtcell_write_bam:
        Read ccs reads from multiple bam files of one SMRTCell,
        select records of a few zmws from each file, and
        write to an output bam file.
        """
        movies = ["m140802_021814_42161_c110036822550000001823106706241500_s1_p0"]

        ccsbamfns = [op.join(self.bigDataDir, "%s.%d.ccs.bam" % (movies[0], n)) for n in (1,2,3)]
        outbamfn = op.join(self.outDir, "test_writebam_ccs_2.bam")

        hns = [58,
               63749,
               113411]
        zmws = ["%s/%s" % (movie, hn) for (movie, hn) in zip(movies, hns)]
        expected_num_ccs_reads = 13821 + 13457 + 17016

        _verify_write_compare_ccs(self, ccsbamfns, zmws, outbamfn,
                                  expected_movies=movies,
                                  expected_len=expected_num_ccs_reads)

    def test_read_ccs_from_multiple_smrtcells_write_one_bam(self):
        """
        test_read_ccs_from_multiple_smrtcells_write_one_bam:
        Read ccs reads from multiple bam files of multiple SMRTCells,
        select records of a few zmws from each movie, and
        write to an output bam file.
        """
        fofn = op.join(self.bigDataDir, "ccsbam.fofn")
        movies = ["m140802_021814_42161_c110036822550000001823106706241500_s1_p0",
                  "m140802_043938_42161_c110036822550000001823106706241501_s1_p0",
                  "m140802_070303_42161_c110036822550000001823106706241502_s1_p0"]

        ccsbamfns = [f for f in readFofn(fofn)]
        outbamfn = op.join(self.outDir, "test_writebam_ccs_3.bam")

        m1, m2, m3 = movies
        hn1 = [54434, 54440, 54493, 80328, 163395]
        hn2 = [1784, 7201, 40789, 79704, 152904]
        hn3 = [705, 4838, 40197, 84197, 126136]
        zmws = ["%s/%d" % (m1, hn) for hn in hn1] + \
               ["%s/%d" % (m2, hn) for hn in hn2] + \
               ["%s/%d" % (m3, hn) for hn in hn3]
        expected_num_ccs_reads = 13821 + 13457 + 17016 + 13705 + 13581 + 14900 + 13477 + 12238 + 11318 #123513

        _verify_write_compare_ccs(self, ccsbamfns, zmws, outbamfn,
                                  expected_movies=movies,
                                  expected_len=expected_num_ccs_reads)


    def test_read_ccs_multiple_movies_one_bam(self):
        """
        Test for sane BamCollection.__getitem__() behavior when a .bam file
        contains multiple read groups.
        """
        dataset_xml = op.join(self.ioDataDir,
                              "ccs_multi_movie.consensusreadset.xml")
        bc = BamCollection(dataset_xml)
        with ConsensusReadSet(dataset_xml) as ds:
            for read in ds:
                self.assertEqual(bc[read.qName].readName, read.qName)
