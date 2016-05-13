"""
Test tool contract with cluster chunking interfaces.
"""
import shutil
import unittest
import logging
import cPickle
import re
import os.path as op
import os

from pbcore.io import ContigSet, FastaReader
from pbcore.util.Process import backticks
from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.pb_io.report import load_report_from_json
import pbcommand.testkit.core

from pbtranscript.ice.IceFiles import IceFiles
from pbtranscript.Utils import mknewdir, mkdir, as_contigset
from pbtranscript.separate_flnc import SeparateFLNCBySize
from pbtranscript.tasks.TPickles import *
from test_setpath import OUT_DIR, DATA_DIR

logging.basicConfig(level=logging.DEBUG)

MNT_DATA = "/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data"
SUBREADS_DATASET = "subreads.subreadset.xml"
CCS_DATASET = "ccs.consensusreadset.xml"
FLNC_DATASET = "isoseq_flnc.contigset.xml"
NFL_DATASET = "isoseq_nfl.contigset.xml"

test_name = "test_tool_contract_chunks"
data_dir = op.join(MNT_DATA, test_name)
subreads_ds = op.join(data_dir, SUBREADS_DATASET)
ccs_ds = op.join(data_dir, CCS_DATASET)
flnc_ds = op.join(data_dir, FLNC_DATASET)
nfl_ds = op.join(data_dir, NFL_DATASET)

done_txt = op.join(data_dir, "done.txt")
cluster_chunks_pickle = op.join(data_dir, "cluster_chunks.pickle")
partial_chunks_pickle = op.join(data_dir, "partial_chunks.pickle")
polish_chunks_pickle = op.join(data_dir, "polish_chunks.pickle")

BIN_NAMES = ["0to1kb_part0", "1to2kb_part0", "2to3kb_part0",
             "3to4kb_part0", "4to5kb_part0"]

CLUSTER_OUT_DIRS = [op.join(data_dir, bin_name, "cluster_out")
                    for bin_name in BIN_NAMES]
N_NFL_CHUNKS = 2
HQ_ISOFORMS_FNS = ["all_quivered_hq.100_30_0.99.fasta",
                   "all_quivered_hq.100_30_0.99.fastq"]
LQ_ISOFORMS_FNS = ["all_quivered_lq.fasta",
                   "all_quivered_lq.fastq"]


def call_separate_flnc(flnc_ds, out_dir, out_pickle):
    """Separate FLNC reads in flnc_ds, and save output to out_dir.
    e.g., out_dir/0to1kb_part0/, out_dir/1to2kb_part0/
    """
    with SeparateFLNCBySize(flnc_filename=flnc_ds, root_dir=out_dir,
                            out_pickle=out_pickle) as obj:
        obj.run()


def make_pickle(in_pickle, out_pickle, root_dir,
                copy_consensus_isoforms=False,
                copy_flnc_pickle=False,
                copy_nfl_pickle=False,
                copy_quivered=False):
    """
    Copy cluster_out_dir in in_pickle to {root_dir}/bin_name/cluster_out/
    """
    mkdir(root_dir)
    def make_flnc(in_flnc, root_dir):
        bin_name = op.basename(op.dirname(in_flnc))
        flnc_name = op.basename(in_flnc)

        assert in_flnc.endswith(".contigset.xml")
        in_flnc_fa = in_flnc.replace(".contigset.xml", ".fasta")
        new_flnc = op.join(root_dir, bin_name, flnc_name)
        new_flnc_fa = new_flnc.replace(".contigset.xml", ".fasta")

        print "new_flnc = %s" % new_flnc
        shutil.copy(in_flnc_fa, new_flnc_fa)
        as_contigset(new_flnc_fa, new_flnc)

    def make_cluster_out_dir(in_dir, root_dir):
        bin_name = op.basename(op.dirname(in_dir))
        new_dir = op.join(root_dir, bin_name, "cluster_out") #e.g., root_dir/0to1kb_part0/cluster_out
        mkdir(new_dir)
        return new_dir

    def _cp(task, new_task, copied_files, copy_consensus_isoforms=copy_consensus_isoforms,
            copy_flnc_pickle=copy_flnc_pickle, copy_nfl_pickle=copy_nfl_pickle):
        """Copy task.files to new_task.files."""
        if copy_consensus_isoforms is True and new_task.consensus_isoforms_file not in copied_files:
            shutil.copy(task.consensus_isoforms_file, new_task.consensus_isoforms_file)
            copied_files[new_task.consensus_isoforms_file] = True
        if copy_flnc_pickle is True and new_task.flnc_pickle not in copied_files:
            mkdir(op.dirname(new_task.flnc_pickle))
            shutil.copy(task.flnc_pickle, new_task.flnc_pickle)
            copied_files[new_task.flnc_pickle] = True
        if copy_nfl_pickle is True and new_task.nfl_pickle not in copied_files:
            mkdir(op.dirname(new_task.nfl_pickle))
            shutil.copy(task.nfl_pickle, new_task.nfl_pickle)
            copied_files[new_task.nfl_pickle] = True


    print "making pickle from in_pickle %s to out_pickle %s, root_dir %s" % \
            (in_pickle, out_pickle, root_dir)

    p = ChunkTasksPickle.read(in_pickle)
    assert len(p) > 0
    if all([isinstance(task, ClusterChunkTask) for task in p]):
        outp = ChunkTasksPickle()
        copied_files = dict()
        for task in p:
            cluster_out_dir = make_cluster_out_dir(task.cluster_out_dir, root_dir)
            print "new_cluster_out_dir is %s" % cluster_out_dir
            #flnc_file = make_flnc(task.flnc_file)
            new_task = ClusterChunkTask(task.cluster_bin_index,
                                        task.flnc_file,
                                        cluster_out_dir)
            _cp(task=task, new_task=new_task, copied_files=copied_files,
                copy_consensus_isoforms=copy_consensus_isoforms,
                copy_flnc_pickle=copy_flnc_pickle, copy_nfl_pickle=copy_nfl_pickle)
            outp.append(new_task)
        outp.write(out_pickle)
    elif all([isinstance(task, PartialChunkTask) for task in p]):
        outp = ChunkTasksPickle()
        copied_files = dict()
        for task in p:
            cluster_out_dir = make_cluster_out_dir(task.cluster_out_dir, root_dir)
            print "new_cluster_out_dir is %s" % cluster_out_dir
            #flnc_file = make_flnc(task.flnc_file)
            new_task = PartialChunkTask(task.cluster_bin_index,
                                        task.flnc_file,
                                        cluster_out_dir,
                                        nfl_file=task.nfl_file,
                                        nfl_index=task.nfl_index,
                                        n_nfl_chunks=task.n_nfl_chunks)
            _cp(task=task, new_task=new_task, copied_files=copied_files,
                copy_consensus_isoforms=copy_consensus_isoforms,
                copy_flnc_pickle=copy_flnc_pickle, copy_nfl_pickle=copy_nfl_pickle)
            outp.append(new_task)
        outp.write(out_pickle)
    elif all([isinstance(task, PolishChunkTask) for task in p]):
        outp = ChunkTasksPickle()
        copied_files = dict()
        for task in p:
            cluster_out_dir = make_cluster_out_dir(task.cluster_out_dir, root_dir)
            print "new_cluster_out_dir is %s" % cluster_out_dir
            #flnc_file = make_flnc(task.flnc_file)
            new_task = PolishChunkTask(task.cluster_bin_index,
                                       task.flnc_file,
                                       cluster_out_dir,
                                       polish_index=task.polish_index,
                                       n_polish_chunks=task.n_polish_chunks)
            mkdir(op.dirname(new_task.nfl_pickle))
            # always copy nfl_pickle for PolishChunkTask
            assert copy_nfl_pickle is True

            _cp(task=task, new_task=new_task, copied_files=copied_files,
                copy_consensus_isoforms=copy_consensus_isoforms,
                copy_flnc_pickle=copy_flnc_pickle, copy_nfl_pickle=copy_nfl_pickle)
            dst_dir =op.join(cluster_out_dir, "quivered")
            if copy_quivered is True and dst_dir not in copied_files:
                if op.exists(dst_dir):
                    shutil.rmtree(dst_dir)
                shutil.copytree(op.join(task.cluster_out_dir, "quivered"), dst_dir)
                copied_files[dst_dir] = True
            outp.append(new_task)
        outp.write(out_pickle)
    else:
        assert False


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestCreateChunks(pbcommand.testkit.PbTestApp):
    """Call python -m pbtranscript.tasks.create_chunks --resolved-tool-contract rtc.json"""
    out_dir = op.join(OUT_DIR, "test_create_chunks")
    mknewdir(out_dir)
    separate_flnc_pickle = op.join(out_dir, "separate_flnc.pickle")
    call_separate_flnc(flnc_ds=flnc_ds, out_dir=out_dir,
                       out_pickle=separate_flnc_pickle)

    DRIVER_BASE = "python -m pbtranscript.tasks.create_chunks"
    INPUT_FILES = [separate_flnc_pickle,  # input 0, separate_flnc.pickle
                   nfl_ds]  # input 1, nfl.xml

    def run_after(self, rtc, output_dir):
        print rtc.task.output_files[0]
        print rtc.task.output_files[1]
        print rtc.task.output_files[2]

@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestClusterBins(pbcommand.testkit.PbTestApp):
    """Call python -m pbtranscript.tasks.cluster_bins --resolved-tool-contract rtc.json"""
    out_dir = op.join(OUT_DIR, "test_cluster_bins")
    mknewdir(out_dir)

    out_cluster_chunks_pickle = op.join(out_dir, "cluster_chunks.pickle")
    make_pickle(in_pickle=cluster_chunks_pickle,
                out_pickle=out_cluster_chunks_pickle,
                root_dir=out_dir)

    DRIVER_BASE = "python -m pbtranscript.tasks.cluster_bins"
    INPUT_FILES = [out_cluster_chunks_pickle,  # input 0, cluster_chunks.pickle
                   ccs_ds] # idx 1, ccs

    def run_after(self, rtc, output_dir):
        out_dir = op.join(OUT_DIR, "test_cluster_bins")
        cluster_out_dirs = [op.join(out_dir, bin_name, "cluster_out")
                            for bin_name in BIN_NAMES]
        self.assertTrue(op.exists(rtc.task.output_files[0]))
        out_consensus_isoforms = [op.join(d, "output", "final.consensus.fasta") for d in cluster_out_dirs]
        print out_consensus_isoforms
        self.assertTrue(all([op.exists(f) for f in out_consensus_isoforms]))


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestIcePartialClusterBins(pbcommand.testkit.PbTestApp):
    """Call python -m pbtranscript.tasks.ice_partial_cluster_bins --resolved-tool-contract rtc.json"""
    out_dir = op.join(OUT_DIR, "test_ice_partial_cluster_bins")
    mknewdir(out_dir)

    out_partial_chunks_pickle = op.join(out_dir, "partial_chunks.pickle")
    make_pickle(in_pickle=partial_chunks_pickle,
                out_pickle=out_partial_chunks_pickle,
                root_dir=out_dir, copy_consensus_isoforms=True)

    DRIVER_BASE = "python -m pbtranscript.tasks.ice_partial_cluster_bins"
    INPUT_FILES = [out_partial_chunks_pickle,  # input 0, partial_chunk.pickle
                   done_txt, # idx 1, sentinel file
                   ccs_ds] # idx 2, ccs

    def run_after(self, rtc, output_dir):
        self.assertTrue(op.exists(rtc.task.output_files[0]))

        out_dir = op.join(OUT_DIR, "test_ice_partial_cluster_bins")
        cluster_out_dirs = [op.join(out_dir, bin_name, "cluster_out")
                            for bin_name in BIN_NAMES]
        out_pickles = [IceFiles(prog_name="", root_dir=d).nfl_pickle_i(i=i)
                       for d in cluster_out_dirs for i in range(N_NFL_CHUNKS)]
        print "output scattered nfl pickles are %s" % out_pickles
        self.assertTrue(all([op.exists(f) for f in out_pickles]))


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestGatherIcePartialPickle(pbcommand.testkit.PbTestApp):
    """Call python -m pbtranscript.tasks.gather_ice_partial_cluster_bins_pickle --resolved-tool-contract rtc.json"""
    out_dir = op.join(OUT_DIR, "test_gather_ice_partial_cluster_bins_pickle")
    mknewdir(out_dir)

    out_partial_chunks_pickle = op.join(out_dir, "partial_chunks.pickle")
    make_pickle(in_pickle=partial_chunks_pickle,
                out_pickle=out_partial_chunks_pickle,
                root_dir=out_dir,
                copy_consensus_isoforms=True,
                copy_nfl_pickle=True)

    DRIVER_BASE = "python -m pbtranscript.tasks.gather_ice_partial_cluster_bins_pickle"
    INPUT_FILES = [out_partial_chunks_pickle,  # input 0, partial_chunk.pickle
                   done_txt] # idx 1, sentinel file

    def run_after(self, rtc, output_dir):
        self.assertTrue(op.exists(rtc.task.output_files[0]))

        out_dir = op.join(OUT_DIR, "test_gather_ice_partial_cluster_bins_pickle")
        cluster_out_dirs = [op.join(out_dir, bin_name, "cluster_out")
                            for bin_name in BIN_NAMES]
        out_pickles = [IceFiles(prog_name="", root_dir=d).nfl_all_pickle_fn
                       for d in cluster_out_dirs]
        print "output nfl pickles are %s" % out_pickles
        self.assertTrue(all([op.exists(f) for f in out_pickles]))


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestIcePolishClusterBins(pbcommand.testkit.PbTestApp):
    """Call python -m pbtranscript.tasks.ice_polish_cluster_bins --resolved-tool-contract rtc.json"""
    out_dir = op.join(OUT_DIR, "test_ice_polish_cluster_bins")
    mknewdir(out_dir)

    out_polish_chunks_pickle = op.join(out_dir, "polish_chunks.pickle")
    make_pickle(in_pickle=polish_chunks_pickle,
                out_pickle=out_polish_chunks_pickle,
                root_dir=out_dir,
                copy_consensus_isoforms=True,
                copy_flnc_pickle=True,
                copy_nfl_pickle=True)

    DRIVER_BASE = "python -m pbtranscript.tasks.ice_polish_cluster_bins"
    INPUT_FILES = [out_polish_chunks_pickle,  # input 0, polish_chunk.pickle
                   done_txt,  # idx 1, sentinel file
                   subreads_ds ] # idx 2, subreads.bam

    def run_after(self, rtc, output_dir):
        self.assertTrue(op.exists(rtc.task.output_files[0]))

@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestGatherPolishedIsoforms(pbcommand.testkit.PbTestApp):
    """Call python -m pbtranscript.tasks.gather_polished_isoforms_in_each_bin --resolved-tool-contract rtc.json"""
    out_dir = op.join(OUT_DIR, "test_gather_polished_isoforms_in_each_bin")
    mknewdir(out_dir)

    out_polish_chunks_pickle = op.join(out_dir, "polish_chunks.pickle")
    make_pickle(in_pickle=polish_chunks_pickle,
                out_pickle=out_polish_chunks_pickle,
                root_dir=out_dir,
                copy_consensus_isoforms=True,
                copy_flnc_pickle=True,
                copy_nfl_pickle=True,
                copy_quivered=True)

    DRIVER_BASE = "python -m pbtranscript.tasks.gather_polished_isoforms_in_each_bin"
    INPUT_FILES = [out_polish_chunks_pickle,  # input 0, polish_chunk.pickle
                   done_txt]  # idx 1, sentinel file

    def run_after(self, rtc, output_dir):
        self.assertTrue(op.exists(rtc.task.output_files[0]))

        out_dir = op.join(OUT_DIR, "test_gather_polished_isoforms_in_each_bin")
        cluster_out_dirs = [op.join(out_dir, bin_name, "cluster_out")
                            for bin_name in BIN_NAMES]
        out_hq_fns = [op.join(d, fn)
                      for d in cluster_out_dirs for fn in HQ_ISOFORMS_FNS]
        print "out_hq_fns %s" % out_hq_fns
        self.assertTrue(all([op.exists(f) for f in out_hq_fns]))

        out_lq_fns = [op.join(d, fn)
                       for d in cluster_out_dirs for fn in LQ_ISOFORMS_FNS]
        print "out_lq_fns %s" % out_lq_fns
        self.assertTrue(all([op.exists(f) for f in out_lq_fns]))

        print "out_lq_fa %s is not empty" % out_lq_fns[0]
        n = len([r for r in FastaReader(out_lq_fns[0])])
        self.assertTrue(n > 0)

        out_logs = [IceFiles(prog_name="", root_dir=d).submitted_quiver_jobs_log
                    for d in cluster_out_dirs]
        print "out_logs %s" % out_logs
        self.assertTrue(all([op.exists(f) for f in out_logs]))


@unittest.skipUnless(op.isdir(MNT_DATA), "Missing %s" % MNT_DATA)
class TestCombineClusterBins(pbcommand.testkit.PbTestApp):
    """Call python -m pbtranscript.tasks.combine_cluster_bins --resolved-tool-contract rtc.json"""
    out_dir = op.join(OUT_DIR, "test_combine_cluster_bins")
    mknewdir(out_dir)

    out_cluster_chunks_pickle = op.join(out_dir, "cluster_chunks.pickle")
    make_pickle(in_pickle=cluster_chunks_pickle,
                out_pickle=out_cluster_chunks_pickle,
                root_dir=out_dir,
                copy_consensus_isoforms=True,
                copy_flnc_pickle=True,
                copy_nfl_pickle=True)

    cluster_out_dirs = [op.join(out_dir, bin_name, "cluster_out")
                        for bin_name in BIN_NAMES]

    for D, d in zip(CLUSTER_OUT_DIRS, cluster_out_dirs):
        polish_log = op.join("log", "submitted_quiver_jobs.txt")
        shutil.copy(op.join(D, polish_log), op.join(d, polish_log))
        for fn in HQ_ISOFORMS_FNS + LQ_ISOFORMS_FNS:
            shutil.copy(op.join(D, fn), op.join(d, fn))

    DRIVER_BASE = "python -m pbtranscript.tasks.combine_cluster_bins"
    INPUT_FILES = [out_cluster_chunks_pickle,  # input 0, cluster_chunk.pickle
                   done_txt]  # idx 1, sentinel file

    def run_after(self, rtc, output_dir):
        self.assertTrue(op.exists(rtc.task.output_files[i]) for i in range(7))

        out_dir = op.join(OUT_DIR, "test_combine_cluster_bins")
        cluster_out_dirs = [op.join(out_dir, bin_name, "cluster_out")
                            for bin_name in BIN_NAMES]

        combined_lq_cs = rtc.task.output_files[5]
        print "combined_lq_fa %s must not be empty" % combined_lq_cs
        n = len([r for r in ContigSet(combined_lq_cs)])
        self.assertTrue(n > 0)

        out_logs = [IceFiles(prog_name="", root_dir=d).submitted_quiver_jobs_log
                    for d in cluster_out_dirs]
        print "out_logs %s" % out_logs
        self.assertTrue(all([op.exists(f) for f in out_logs]))
