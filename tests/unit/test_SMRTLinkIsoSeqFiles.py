#!/usr/bin/env python

from nose.tools import assert_equal
from pbtranscript.io import SMRTLinkIsoSeqFiles

def test_SMRTLinkIsoSeqFiles():
    f = SMRTLinkIsoSeqFiles('/root')
    assert_equal(f.tasks_dir, '/root/tasks')
    assert_equal(f.combined_dir, '/root/tasks/pbtranscript.tasks.separate_flnc-0/combined')
    assert_equal(f.ccs_ds, '/root/tasks/pbcoretools.tasks.gather_ccsset-1/file.consensusreadset.xml')
    assert_equal(f.ccs_fa_gz, '/root/tasks/pbcoretools.tasks.bam2fasta_ccs-0/ccs.gz')
    assert_equal(f.flnc_gather_dir, '/root/tasks/pbcoretools.tasks.gather_contigset-2')
    assert_equal(f.isoseq_flnc_ds, '/root/tasks/pbcoretools.tasks.gather_contigset-2/file.contigset.xml')
    assert_equal(f.isoseq_flnc_fa, '/root/tasks/pbcoretools.tasks.gather_contigset-2/file.contigset.fasta')

    assert_equal(f.nfl_gather_dir, '/root/tasks/pbcoretools.tasks.gather_contigset-3')
    assert_equal(f.isoseq_nfl_ds, '/root/tasks/pbcoretools.tasks.gather_contigset-3/file.contigset.xml')
    assert_equal(f.isoseq_nfl_fa, '/root/tasks/pbcoretools.tasks.gather_contigset-3/file.contigset.fasta')
    assert_equal(f.draft_gather_dir, '/root/tasks/pbcoretools.tasks.gather_contigset-1')
    assert_equal(f.isoseq_draft_ds, '/root/tasks/pbcoretools.tasks.gather_contigset-1/file.contigset.xml')
    assert_equal(f.isoseq_draft_fa, '/root/tasks/pbcoretools.tasks.gather_contigset-1/file.contigset.fasta')
    assert_equal(f.hq_isoforms_fa, '/root/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.polished_hq.fasta')
    assert_equal(f.hq_isoforms_fq, '/root/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.polished_hq.fastq')
    assert_equal(f.lq_isoforms_fa, '/root/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.polished_lq.fasta')
    assert_equal(f.lq_isoforms_fq, '/root/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.polished_lq.fastq')
    assert_equal(f.consensus_isoforms_fa, '/root/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.consensus_isoforms.fasta')
    assert_equal(f.ccs_report_json, '/root/tasks/pbreports.tasks.ccs_report-0/ccs_report.json')
    assert_equal(f.classify_report_json, '/root/tasks/pbreports.tasks.isoseq_classify-0/isoseq_classify_report.json')
    assert_equal(f.cluster_report_csv, '/root/tasks/pbtranscript.tasks.separate_flnc-0/combined/cluster_report.csv')
    assert_equal(f.hq_lq_prefix_pickle, '/root/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.hq_lq_pre_dict.pickle')
