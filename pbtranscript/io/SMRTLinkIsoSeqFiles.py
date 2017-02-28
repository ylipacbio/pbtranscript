#!/usr/bin/env python

"""Access data files in SMRTLink IsoSeq job"""

import logging
import os.path as op
from pbtranscript.Utils import realpath

log = logging.getLogger(__name__)

__all__ = ["SMRTLinkIsoSeqFiles"]

class SMRTLinkIsoSeqFiles(object):

    """Files of SMRTLink IsoSeq job"""

    def __init__(self, root_dir):
        self.root_dir = realpath(root_dir)

    @property
    def tasks_dir(self):
        """Return ${root_dir}/tasks"""
        return op.join(self.root_dir, "tasks")
    @property
    def subreads_xml(self):
        """return subreads eid_subread"""
        def get_subreads_xml_from_jobstdout(job_stdout):
            for l in open(job_stdout, 'r'):
                if "\'eid_subread': " in l:
                    return l.split(':')[1].translate(None, '}\' ').strip()
        try:
            return get_subreads_xml_from_jobstdout(op.join(self.root_dir, 'job.stdout'))
        except Exception as e:
            raise ValueError("Could not find subreads.xml for smrtlink job %s %s" % (self.root_dir, str(e)))

    @property
    def combined_dir(self):
        """combined dir for isoseq"""
        return op.join(self.tasks_dir, "pbtranscript.tasks.separate_flnc-0", "combined")

    @property
    def ccs_ds(self):
        """file path to ccs dataset xml"""
        return op.join(self.tasks_dir, "pbcoretools.tasks.gather_ccsset-1/file.consensusreadset.xml")

    @property
    def ccs_fa_gz(self):
        """file path to ccs fasta gz"""
        return op.join(self.tasks_dir, "pbcoretools.tasks.bam2fasta_ccs-0/ccs.gz")

    @property
    def isoseq_flnc_ds(self):
        """file path to isoseq_flnc contigset xml"""
        return op.join(self.tasks_dir, "pbcoretools.tasks.gather_contigset-3/file.contigset.xml")

    @property
    def isoseq_nfl_ds(self):
        """file path to isoseq_nfl contigset xml"""
        return op.join(self.tasks_dir, "pbcoretools.tasks.gather_contigset-2/file.contigset.xml")

    @property
    def isoseq_draft_ds(self):
        """file path to isoseq_draft contigset xml"""
        return op.join(self.tasks_dir, "pbcoretools.tasks.gather_contigset-1/file.contigset.xml")

    @property
    def hq_isoforms_fa(self):
        """file path to hq isoforms.fasta"""
        return op.join(self.combined_dir, "all.polished_hq.fasta")

    @property
    def hq_isoforms_fq(self):
        """file path to hq isoforms.fasta"""
        return op.join(self.combined_dir, "all.polished_hq.fastq")

    @property
    def lq_isoforms_fa(self):
        """file path to lq isoforms.fasta"""
        return op.join(self.combined_dir, "all.polished_lq.fasta")

    @property
    def lq_isoforms_fq(self):
        """file path to lq isoforms.fasta"""
        return op.join(self.combined_dir, "all.polished_lq.fastq")

    @property
    def consensus_isoforms_fa(self):
        """file path to consensus isoforms.fasta"""
        return op.join(self.combined_dir, "all.consensus_isoforms.fasta")

    @property
    def ccs_report_json(self):
        """file path to ccs report.json"""
        return op.join(self.tasks_dir, "pbreports.tasks.ccs_report-0", "ccs_report.json")

    @property
    def classify_report_json(self):
        """file path to classify report.json"""
        return op.join(self.tasks_dir, "pbreports.tasks.isoseq_classify-0", "isoseq_classify_report.json")

    @property
    def cluster_report_json(self):
        """file path to classify report.json"""
        return op.join(self.tasks_dir, "pbreports.tasks.isoseq_cluster-0", "isoseq_cluster_report.json")

    @property
    def cluster_report_csv(self):
        """file path to classify report.csv"""
        return op.join(self.combined_dir, "cluster_report.csv")

    @property
    def hq_lq_prefix_pickle(self):
        """file path to classify report.csv"""
        return op.join(self.combined_dir, "all.hq_lq_pre_dict.pickle")
