import unittest
import os.path as op
import filecmp
import numpy as np
from pbcore.util.Process import backticks
from pbtranscript.ClusterOptions import SgeOptions
from pbtranscript.Utils import make_pbi, mknewdir
from pbtranscript.ice import IceUtils
from pbtranscript.io import *
from pbtranscript.io.PbiBamIO import BamZmwRead
from pbtranscript.ice.ProbModel import ProbFromModel
from pbtranscript.ice.IceUtils import *
from pbtranscript.ice_daligner import DalignerRunner
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR, SIV_STD_DIR

def _check_trimmed_BamZmwRead(trimmed_read, original_read, trim_len):
    """Compare trimmed BamZmwRead with the original Zmw read."""
    assert(isinstance(trimmed_read, BamZmwRead) and
           isinstance(original_read, BamZmwRead))
    qvnames = ["deletionqv", "insertionqv", "mergeqv",
               "substitutionqv", "deletiontag", "ipd"]
    return (trimmed_read.movieName == original_read.movieName and
            trimmed_read.holeNumber == original_read.holeNumber and
            trimmed_read.readStart == original_read.readStart + trim_len and
            trimmed_read.readEnd == original_read.readEnd - trim_len and
            trimmed_read.basecalls() == original_read.basecalls()[trim_len:(len(original_read.basecalls())-trim_len)] and
            np.all([trimmed_read.qv(qvname) == original_read.qv(qvname)[trim_len:(len(original_read.qv(qvname))-trim_len)] for qvname in qvnames]))


def _trim_subreads_and_write(in_bam, zmws, out_bam, trim_len, min_len):
    """Test IceUtils.trim_subreads_and_write."""
    assert(op.exists(in_bam))

    bam_reader = BamCollection(in_bam)

    out_m = IceUtils.trim_subreads_and_write(bam_reader,
            zmws, out_bam, trim_len=trim_len,
            min_len=min_len, ignore_keyerror=False, bam=True)

    assert(out_m == set([zmw.split('/')[0] for zmw in zmws]))

    print "Making pbindex for %s" % out_bam
    make_pbi(out_bam)
    print "pbindex ended"

    out_bam_reader = BamCollection(out_bam)

    for zmw in zmws:
        out_srs = out_bam_reader[zmw].subreads
        in_srs = bam_reader[zmw].subreads
        assert(len(out_srs) == len(in_srs))
        for out_sr, in_sr in zip(out_srs, in_srs):
            assert(_check_trimmed_BamZmwRead(out_sr, in_sr, trim_len))


class Test_ICEUtils(unittest.TestCase):
    """Test IceUtils."""
    def setUp(self):
        """Initialize."""
        self.outDir       = OUT_DIR
        self.dataDir      = DATA_DIR
        self.sivDataDir   = SIV_DATA_DIR
        self.sivStdoutDir = SIV_STD_DIR
        self.moreDir = op.join(self.sivDataDir, "test_IceUtils")

    def test_sanity_check_daligner(self):
        """sanity_check_daligner."""
        self.assertTrue(IceUtils.sanity_check_daligner(scriptDir=self.outDir))

    def test_sanity_check_gcon(self):
        """sanity_check_gcon."""
        self.assertTrue(IceUtils.sanity_check_gcon() == IceUtils.gcon_py)

    @unittest.skipUnless(backticks('qstat')[1] == 0, 'sge disabled')
    def test_sanity_check_sge(self):
        """sanity_check_sge."""
        self.assertTrue(IceUtils.sanity_check_sge(SgeOptions(100), self.outDir))

    def test_ice_fa2fq(self):
        """Test ice_fa2fq extracting qvs from h5/bam to fq."""
        in_fa = op.join(self.sivDataDir, "flnc.fasta")
        ccs_fofn = op.join(self.sivDataDir, "ccs.fofn")
        out_fq = op.join(self.outDir, "test_ice_fa2fq.fastq")
        stdout_fq = op.join(self.sivStdoutDir, "test_ice_fa2fq.fastq")
        # Load QVs from bax.h5
        IceUtils.ice_fa2fq(in_fa, ccs_fofn, out_fq)
        self.assertTrue(filecmp.cmp(out_fq, stdout_fq))

        ccs_bam_fofn = op.join(self.sivDataDir, "ccsbam.fofn")
        out_bam_fq = op.join(self.outDir, "test_ice_fa2fq.bam.fastq")
        # Load QVs from bam
        IceUtils.ice_fa2fq(in_fa, ccs_bam_fofn, out_bam_fq)
        self.assertTrue(filecmp.cmp(out_bam_fq, stdout_fq))

    def test_parsed_read_name(self):
        names = ['m/1234/CCS',
                 'm/1234/0_100_CCS',
                 'm/1234/0_100',
                 'm/1234/100_0_CCS',
                 'm/1234/100_0']
        from pbtranscript.ice.IceUtils import _Parsed_Read_Name
        prns = [_Parsed_Read_Name(name) for name in names]
        self.assertTrue(names == [prn.__repr__() for prn in prns])

    def test_is_blank(self):
        """Test is_blank_sam, and is_blank_bam."""
        blank_bam = op.join(self.moreDir, "blank.bam")
        blank_sam = op.join(self.moreDir, "blank.sam")
        non_blank_bam = op.join(self.moreDir, "non_blank.bam")
        non_blank_sam = op.join(self.moreDir, "non_blank.sam")
        from pbtranscript.ice.IceUtils import is_blank_sam, is_blank_bam
        self.assertTrue(is_blank_bam(blank_bam))
        self.assertTrue(is_blank_sam(blank_sam))
        self.assertFalse(is_blank_bam(non_blank_bam))
        self.assertFalse(is_blank_sam(non_blank_sam))

    def test_concat_bam_header(self):
        """Test concat_bam_header"""
        from pbtranscript.ice.IceUtils import concat_bam_header
        fns = [op.join(self.moreDir, "aligned.%d.bam" % i) for i in range(1, 6)]
        out_fn =  op.join(self.outDir, "test_concat_bam_header.sam")
        stdout_fn =  op.join(self.sivStdoutDir, "test_concat_bam_header.sam")
        concat_bam_header(fns, out_fn)
        self.assertTrue(op.exists(out_fn))
        self.assertTrue(filecmp.cmp(out_fn, stdout_fn))

    def test_concat_bam(self):
        """Test concat_bam, unaligned and aligned."""
        # cat aligned bam files with only one RG, one SN
        fns = [op.join(self.moreDir, "%d.bam" % i) for i in range(1,5)]
        out_fn = op.join(self.outDir, "test_concat_bam_1.bam")
        from pbtranscript.ice.IceUtils import concat_bam
        concat_bam(fns, out_fn)
        self.assertTrue(op.exists(out_fn))

        # cat aligned bam files to a big bam
        fns = [op.join(self.moreDir, "aligned.%d.bam" % i) for i in range(1, 6)]
        out_fn = op.join(self.outDir, "test_concat_bam_2.bam")
        from pbtranscript.Utils import execute
        concat_bam(fns, out_fn)
        self.assertTrue(op.exists(out_fn))

        # convert big bam to sam and compare with std output
        out_sam = out_fn + ".sam"
        stdout_sam = op.join(self.sivStdoutDir, "test_concat_bam_2.sam")
        cmd="samtools view -h %s -o %s" % (out_fn, out_sam)
        execute(cmd=cmd)
        self.assertTrue(filecmp.cmp(out_sam, stdout_sam))

    def test_trim_subreads_and_write(self):
        """
        Test trim_subreads_and_write(reader, in_zmwids, outfile, trim_len, min_len...)
        """
        m = "m131018_081703_42161_c100585152550000001823088404281404_s1_p0"
        hns = [45, 161, 227, 293, 495, 642, 780, 865, 888]
        zmws = ["%s/%s" % (m, hn) for hn in hns]

        in_bam = op.join(self.sivDataDir, "bam.fofn")

        # No trim subreads flanks, in == out
        trim_len, min_len = 0, 0
        out_bam = op.join(self.outDir,
                          "test_trim_subreads_%s_%s.bam" % (trim_len, min_len))
        _trim_subreads_and_write(in_bam=in_bam,
                                 zmws=zmws,
                                 out_bam=out_bam,
                                 trim_len=trim_len,
                                 min_len=min_len)

        # Trim both ends of subreads by 100 bp
        trim_len, min_len = 100, 0
        _filename = "test_trim_subreads_%s_%s.bam" % (trim_len, min_len)
        out_bam = op.join(self.outDir, _filename)
        std_out_bam = op.join(self.sivStdoutDir, _filename)
        _trim_subreads_and_write(in_bam=in_bam,
                                 zmws=zmws,
                                 out_bam=out_bam,
                                 trim_len=trim_len,
                                 min_len=min_len)

    def _test_daligner_against_ref(self, test_name, use_sge, sge_opts,
                                   prob_model_from="fake"):
        """Test daligner_against_ref with and without using sge."""
        copy_dir = op.join(self.dataDir, "test_daligner_against_ref")
        output_dir = op.join(self.outDir, test_name)
        mknewdir(output_dir)

        qname, tname = "test_daligner_query.fasta", "test_daligner_target.fasta"
        query_filename = op.join(output_dir, qname)
        target_filename = op.join(output_dir, tname)

        prob_model = None
        if prob_model_from == "fake":
            prob_model = ProbFromModel(0.01, 0.07, 0.06)
        elif prob_model_from == "fastq":
            fastq_fn = op.join(copy_dir, "test_daligner_reads.fastq")
            prob_model = ProbFromFastq(fastq_fn)
        else:
            self.assertTrue(False)

        qver_get_func = prob_model.get_smoothed
        qvmean_get_func = prob_model.get_mean

        dummy_o, c, dummy_m = backticks("cp %s %s" % (op.join(copy_dir, qname), query_filename))
        self.assertTrue(c == 0)

        dummy_o, c, dummy_m = backticks("cp %s %s" % (op.join(copy_dir, tname), target_filename))
        self.assertTrue(c == 0)

        old_dir = os.getcwd()
        os.chdir(output_dir)

        runner = DalignerRunner(query_filename=query_filename,
                                target_filename=target_filename,
                                is_FL=True, same_strand_only=True,
                                use_sge=use_sge, sge_opts=sge_opts)
        runner.run(output_dir=op.join(self.outDir, test_name))

        hits = []

        for la4ice_filename in runner.la4ice_filenames:
            hits.extend(daligner_against_ref(query_dazz_handler=runner.query_dazz_handler,
                                             target_dazz_handler=runner.target_dazz_handler,
                                             la4ice_filename=la4ice_filename,
                                             is_FL=True, sID_starts_with_c=False,
                                             qver_get_func=qver_get_func,
                                             qvmean_get_func=qvmean_get_func))
        # Num of hits may change when daligner or parameters change.
        self.assertTrue(len(hits), 706)
        self.assertEqual(str(hits[0]),
                         "m54007_160109_025449/27984844/29_646_CCS/0_617 aligns to m54007_160109_025449/28836279/631_54_CCS")
        os.chdir(output_dir)

    def test_daligner_against_ref(self):
        """Test daligner_against_ref() using fake prob model."""
        test_name = "test_daligner_against_ref"
        self._test_daligner_against_ref(test_name=test_name, use_sge=False, sge_opts=None)

    def test_daligner_against_ref_prob_model_fastq(self):
        """Test daligner_against_ref() using prob model from FASTQ"""
        test_name = "test_daligner_against_ref_prob_model_fastq"
        self._test_daligner_against_ref(test_name=test_name, use_sge=False, sge_opts=None,
                                        prob_model_from="fastq")

    @unittest.skipUnless(backticks('qstat')[1] == 0, "sge disabled")
    def test_daligner_against_ref_use_sge(self):
        """Test daligner_against_ref() using fake prob model on sge."""
        test_name = "test_daligner_against_ref_use_sge"
        self._test_daligner_against_ref(test_name=test_name, use_sge=True, sge_opts=SgeOptions(100))

    def test_num_reads_in_fasta(self):
        """Test num_reads_in_fasta"""
        in_fa = op.join(self.sivDataDir, "flnc.fasta")
        in_xml = op.join(self.sivDataDir, "test_tool_contract_chunks/isoseq_flnc.contigset.xml")
        self.assertEqual(226, num_reads_in_fasta(in_fa))
        self.assertEqual(161, num_reads_in_fasta(in_xml))
