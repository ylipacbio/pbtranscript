
"""
Test tool contract interfaces.  This also doubles as a test for dataset support
for both .bam and .fasta inputs.
"""

import unittest
import logging
import cPickle
import re
import os.path as op
import os


from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.pb_io.report import load_report_from_json
import pbcommand.testkit.core
from pbcore.io import ContigSet

TEST_DIR = op.dirname(op.dirname(__file__))
ROOT_DIR = op.dirname(TEST_DIR)
DATA = op.join(TEST_DIR, "data")
MNT_DATA = "/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data"
SUBREADS_DATASET = "m131018_081703_42161_c100585152550000001823088404281404_s1_p0.1.subreadset.xml"
CCS_DATASET = "m131018_081703_42161_c100585152550000001823088404281404_s1_p0.1.consensusreadset.xml"
FLNC_DATASET = "isoseq_flnc.contigset.xml"
NFL_DATASET = "isoseq_nfl.contigset.xml"

logging.basicConfig(level=logging.DEBUG)


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestPbtranscriptClassify(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m pbtranscript.tasks.classify"
    REQUIRES_PBCORE = True
    MAX_NPROC = 8
    RESOLVED_NPROC = 8
    INPUT_FILES = [op.join(MNT_DATA, "bam", CCS_DATASET)]

    def run_after(self, rtc, output_dir):
        rpt = None
        uuids = []
        for file_name in rtc.task.output_files:
            if file_name.endswith(".json"):
                rpt = load_report_from_json(file_name)
            elif file_name.endswith(".xml"):
                uuids.append(ContigSet(file_name, strict=True).uuid)
            else:
                assert file_name.endswith(".csv")
        self.assertEqual(sorted(rpt._dataset_uuids), sorted(uuids))


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestPbtranscriptCluster(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m pbtranscript.tasks.cluster"
    REQUIRES_PBCORE = True
    MAX_NPROC = 8
    RESOLVED_NPROC = 8
    INPUT_FILES = [
        op.join(MNT_DATA, "bam", FLNC_DATASET),
        op.join(MNT_DATA, "bam", NFL_DATASET),
        op.join(MNT_DATA, "bam", CCS_DATASET),
        op.join(MNT_DATA, "bam", SUBREADS_DATASET),
    ]


class TestPbtranscriptSubset(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m pbtranscript.tasks.subset"
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        op.join(DATA, "test_subset.fasta"),
    ]


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestIcePartial(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m pbtranscript.tasks.ice_partial"
    INPUT_FILES = [
        op.join(MNT_DATA, "sa3", "nfl.contigset.xml"),
        op.join(MNT_DATA, "sa3", "final_consensus.contigset.xml"),
        op.join(MNT_DATA, "sa3", "ccsbam.consensusreadset.xml")
    ]


# FIXME this is too tied into directory structure to be unit-testable
# class TestIceQuiver(pbcommand.testkit.PbTestApp):
#    DRIVER_BASE = "python -m pbtranscript.tasks.ice_quiver"
#    INPUT_FILES = [
#        op.join(MNT_DATA, "sa3", "bam.subreadset.xml"),
#        op.join(MNT_DATA, "sa3", "final_consensus.contigset.xml"),
#        op.join(MNT_DATA, "sa3", "clusters.pickle"),
#        op.join(MNT_DATA, "sa3", "nfl.all.partial_uc.pickle"),
#    ]

# TODO ice_quiver_postprocess


#@unittest.skip("GMAP disabled")
#class TestGMAP(pbcommand.testkit.core.PbTestApp):
#    DRIVER_BASE = "python -m pbtranscript.tasks.gmap"
#    INPUT_FILES = [
#        op.join(MNT_DATA, "sa3", "nfl.contigset.xml"),
#        "/pbi/dept/secondary/siv/references/rat_UCSC/rat_UCSC.referenceset.xml",
#    ]
#    MAX_NPROC = 4
#    RESOLVED_NPROC = 4
#
#    @classmethod
#    def setUpClass(cls):
#        super(TestGMAP, cls).setUpClass()
#        # XXX workaround for hardcoded paths in gmap_build
#        if not "_SMRT_GMAP_BIN" in os.environ:
#            top_dir = op.dirname(op.dirname(op.dirname(ROOT_DIR)))
#            prebuilt_dir = op.join(top_dir, "prebuilt.out")
#            if op.exists(prebuilt_dir):
#                gmap_dir = op.join(prebuilt_dir,
#                                   "gmap/gmap-2014-12-21/ubuntu-1404/gmap_home/bin")
#                os.environ["_SMRT_GMAP_BIN"] = gmap_dir
#                os.environ["PATH"] = ":".join([os.environ["PATH"], gmap_dir])
#            else:
#                raise OSError("Need to define _SMRT_GMAP_BIN")


########################################################################
# CHUNKING SUPPORT
########################################################################

class ContigSetScatterBase(object):
    def run_after(self, rtc, output_dir):
        json_file = rtc.task.output_files[0]
        chunks = load_pipeline_chunks_from_json(json_file)
        n_rec = 0
        with ContigSet(self.INPUT_FILES[0]) as f:
            n_rec = len(f)
        n_rec_chunked = 0
        for chunk in chunks:
            d = chunk.chunk_d
            with ContigSet(d['$chunk.contigset_id']) as cs:
                n_rec_chunked += len([r for r in cs])
            self._check_unchunked_files(d)
        self.assertEqual(n_rec_chunked, n_rec)

    def _check_unchunked_files(self, d):
        return True


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing {d}".format(d=MNT_DATA))
class TestScatterContigSet(ContigSetScatterBase,
                           pbcommand.testkit.core.PbTestScatterApp):
    DRIVER_BASE = "python -m pbtranscript.tasks.scatter_contigset"
    INPUT_FILES = [  # identical to TestIcePartial inputs
        op.join(MNT_DATA, "sa3", "nfl.contigset.xml"),
        op.join(MNT_DATA, "sa3", "flnc.contigset.xml"),
        op.join(MNT_DATA, "sa3", "ccsbam.consensusreadset.xml"),
    ]

    def _check_unchunked_files(self, d):
        self.assertEqual(d['$chunk.ref_contigset_id'], self.INPUT_FILES[1])
        self.assertEqual(d['$chunk.ccsset_id'], self.INPUT_FILES[2])


#@unittest.skip("GMAP disabled")
#class TestScatterContigSetGMAP(ContigSetScatterBase,
#                               pbcommand.testkit.core.PbTestScatterApp):
#    DRIVER_BASE = "python -m pbtranscript.tasks.scatter_contigset_gmap"
#    INPUT_FILES = [  # identical to TestIcePartial inputs
#        op.join(MNT_DATA, "sa3", "nfl.contigset.xml"),
#        "/pbi/dept/secondary/siv/references/rat_UCSC/rat_UCSC.referenceset.xml",
#    ]
#
#    def _check_unchunked_files(self, d):
#        self.assertEqual(d['$chunk.reference_id'], self.INPUT_FILES[1])


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestScatterClusters(pbcommand.testkit.core.PbTestScatterApp):
    DRIVER_BASE = "python -m pbtranscript.tasks.scatter_clusters"
    INPUT_FILES = [
        op.join(MNT_DATA, "sa3", "bam.subreadset.xml"),
        op.join(MNT_DATA, "sa3", "flnc.contigset.xml"),
        op.join(MNT_DATA, "sa3", "clusters.pickle"),
        op.join(MNT_DATA, "sa3", "nfl.all.partial_uc.pickle"),
    ]

    def run_after(self, rtc, output_dir):
        json_file = rtc.task.output_files[0]
        chunks = load_pipeline_chunks_from_json(json_file)
        for chunk in chunks:
            d = chunk.chunk_d
            # the cluster pickle (file index 2) is chunked, the rest not
            self.assertNotEqual(d["$chunk.pickle_id"], self.INPUT_FILES[2])
            self.assertEqual(d["$chunk.subreadset_id"], self.INPUT_FILES[0])
            self.assertEqual(d["$chunk.contigset_id"], self.INPUT_FILES[1])
            self.assertEqual(d["$chunk.nfl_pickle_id"], self.INPUT_FILES[3])


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestGatherPickle(pbcommand.testkit.core.PbTestGatherApp):
    DRIVER_BASE = "python -m pbtranscript.tasks.gather_nfl_pickle"
    INPUT_FILES = [
        op.join(MNT_DATA, "sa3", "chunk_misc", "gathered.chunks.json")
    ]

    def run_after(self, rtc, output_dir):
        pkl = rtc.task.output_files[0]
        with open(pkl, 'rb') as f:
            result = cPickle.load(f)
            assert sorted(result.keys()) == ['nohit', 'partial_uc']


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestSeparateFLNC(pbcommand.testkit.PbTestApp):
    """Call python -m pbtranscript.tasks.separate_flnc --resolved-tool-contract rtc.json"""
    DRIVER_BASE = "python -m pbtranscript.tasks.separate_flnc"
    INPUT_FILES = [op.join(MNT_DATA, "bam", FLNC_DATASET)]

    def run_after(self, rtc, output_dir):
        pkl = rtc.task.output_files[0]
        with open(pkl, 'rb') as f:
            result = cPickle.load(f)
            print 'separate_flnc.pickle is %s' % f
            assert len(result.keys()) == 2
            assert 'sorted_keys' in result.keys()
            assert 'key_to_file' in result.keys()
            assert [op.exists(f) for f in result['key_to_file'].values()]


if __name__ == "__main__":
    unittest.main()
