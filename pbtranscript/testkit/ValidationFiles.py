#!/usr/bin/env python

"""Define isoseq validation files."""

import logging
import os.path as op
from collections import defaultdict
from pbcore.io import ContigSet, FastaReader, FastqReader
from pbtranscript.io import ContigSetReaderWrapper
from pbtranscript.Utils import realpath, rmpath, ln, mkdir
from pbtranscript.testkit.SMRTLinkIsoSeqFiles import SMRTLinkIsoSeqFiles
from pbtranscript.testkit.Utils import consolidate_xml, json_to_attr_dict

FORMATTER = op.basename(__file__) + ':%(levelname)s:'+'%(message)s'
logging.basicConfig(level=logging.INFO, format=FORMATTER)
log = logging.getLogger(__name__)

__all__ = ["ValidationFiles", "ValidationRunner"]

def make_readlength_csv(fasta_fn, csv_fn):
    """Make read length csv file for fasta file"""
    log.info("Making read length csv file %s from %s", csv_fn, fasta_fn)
    rmpath(csv_fn)
    with open(csv_fn, 'w') as writer:
        writer.write("'name'\t'readlength'\n")
        for read in ContigSetReaderWrapper(fasta_fn):
            writer.write('%s\t%s\n' % (read.name.split()[0], len(read.sequence)))


def make_polymerase_readlength_csv(subreads_xml, csv_fn):
    """Use position of the last base in subreads to approximate
    polymerase read lengths, and write them to csv"""
    try:
        from pbcore.io import SubreadSet
        # Assuming all reads from one zmw are grouped
        ds = SubreadSet(subreads_xml)
        # dict: movie id --> movie name
        movie_id_to_name = {movie_id:movie_name for movie_name, movie_id in ds.movieIds.iteritems()}

        ends = defaultdict(lambda: 0)
        for movie_id, zmw, end in zip(ds.index['qId'], ds.index['holeNumber'], ds.index['qEnd']):
            ends[(movie_id, zmw)] = max(end,  ends[(movie_id, zmw)])

        writer = open(csv_fn, 'w')
        writer.write("'name'\t'readlength'\n")
        for movie_id, zmw in sorted(ends.keys()):
            writer.write("%s/%s\t%s\n" % (movie_id_to_name[movie_id], zmw, ends[(movie_id, zmw)]))
        writer.close()
    except Exception as e:
        log.warning("Could not obtain polymerase read length! %s", str(e))


def tolist(x):
    """convert x to a list"""
    return x if isinstance(x, list) else [x]


class ValidationFiles(object):

    """Simple class for all isoseq validation output files."""

    def __init__(self, root_dir):
        self.root_dir = realpath(root_dir)

    @property
    def csv_dir(self):
        """dir for keeping all csv output"""
        return op.join(self.root_dir, 'csv')

    @property
    def fasta_dir(self):
        """dir for keeping all isoseq fasta files"""
        return op.join(self.root_dir, 'fasta')

    @property
    def validation_report_csv(self):
        """file path to validation report"""
        return op.join(self.root_dir, "isoseq_rc0_validation_report.csv")

    @property
    def sirv_report_txt(self):
        """Report for SIRV isoforms"""
        return op.join(self.root_dir, 'SIRV_evaluation_summary.txt')

    @property
    def hg_report_txt(self):
        """Report for hg isoforms"""
        return op.join(self.root_dir, 'human_evaluation_summary.txt')

    @property
    def gencode_gtf(self):
        """file path to human gencode_gtf file"""
        return op.join(self.root_dir, "gencode.annotation.gtf")

    @property
    def sirv_truth_dir(self):
        """dir for keeping sirv ground truth files,
        including touse.gff, touse.ground.txt, touse.count.txt"""
        return op.join(self.root_dir, "SIRV")

    @property
    def readme_txt(self):
        """README file."""
        return op.join(self.root_dir, 'README_.txt')

    @property
    def chain_sample_dir(self):
        """chain sample dir"""
        return op.join(self.root_dir, 'chain_sample')

    @property
    def chain_sample_config(self):
        """chain sample config"""
        return op.join(self.chain_sample_dir, 'chain_sample.contig')

    @property
    def chained_ids_txt(self):
        """all_samples.chained_ids.txt"""
        return op.join(self.chain_sample_dir, 'all_samples.chained_ids.txt')

    @property
    def polymerase_readlength_csv(self):
        """Return polymerase read length csv"""
        return op.join(self.csv_dir, "polymerase_readlength.csv")

    @property
    def ccs_readlength_csv(self):
        """file path to CCS read length csv"""
        return op.join(self.csv_dir, "ccs_readlength.csv")

    @property
    def flnc_readlength_csv(self):
        """file path to FLNC read length csv"""
        return op.join(self.csv_dir, "flnc_readlength.csv")

    @property
    def nfl_readlength_csv(self):
        """file path to NFL read length csv"""
        return op.join(self.csv_dir, "nfl_readlength.csv")

    @property
    def consensus_isoforms_readlength_csv(self):
        """file path to consensus isoforms"""
        return op.join(self.csv_dir, "consensus_isoforms_readlength.csv")

    @property
    def hq_readlength_csv(self):
        """file path to HQ read length"""
        return op.join(self.csv_dir, "hq_readlength.csv")

    @property
    def lq_readlength_csv(self):
        """file path to LQ isoforms read length csv"""
        return op.join(self.csv_dir, "lq_readlength.csv")

    @property
    def collapsed_isoforms_readlength_csv(self):
        """file path to collapsed isoforms"""
        return op.join(self.csv_dir, "collapsed_isoforms_readlength.csv")

    @property
    def ccs_fa_gz(self):
        """file path to concatenated ccs.fasta.gz, no need to
        decompress gz because pbcore.io.FastaReader can read gz files."""
        return op.join(self.fasta_dir, "ccs.fasta.gz")

    @property
    def isoseq_flnc_fa(self):
        """file path to concatenated isoseq_flnc.fasta"""
        return op.join(self.fasta_dir, "isoseq_flnc.fasta")

    @property
    def isoseq_nfl_fa(self):
        """file path to concatenated isoseq_flnc.fasta"""
        return op.join(self.fasta_dir, "isoseq_flnc.fasta")

    @property
    def consensus_isoforms_fa(self):
        """file path to consensus_isoforms.fasta"""
        return op.join(self.fasta_dir, "consensus_isoforms.fasta")

    @property
    def hq_isoforms_fa(self):
        """file path to hq isoforms.fasta"""
        return op.join(self.fasta_dir, "hq_isoforms.fasta")

    @property
    def hq_isoforms_fq(self):
        """file path to hq isoforms.fastq"""
        return op.join(self.fasta_dir, "hq_isoforms.fastq")

    @property
    def lq_isoforms_fa(self):
        """file path to lq isoforms.fasta"""
        return op.join(self.fasta_dir, "lq_isoforms.fasta")

    @property
    def lq_isoforms_fq(self):
        """file path to lq isoforms.fastq"""
        return op.join(self.fasta_dir, "lq_isoforms.fastq")

    @property
    def collapse_to_hg_dir(self):
        """dir for saving files collapsing HQ isoforms to human reference"""
        return op.join(self.root_dir, "collapse_to_hg")

    @property
    def collapse_to_sirv_dir(self):
        """dir for saving files collapsing HQ isoforms to SIRV reference"""
        return op.join(self.root_dir, "collapse_to_sirv")

    @property
    def collapsed_to_hg_rep_fq(self):
        """Representative isoforms collapsed to human genome"""
        return op.join(self.collapse_to_hg_dir, "touse.fastq")

    @property
    def collapsed_to_sirv_rep_fq(self):
        """Representative isoforms collapsed to sirv genome"""
        return op.join(self.collapse_to_sirv_dir, "touse.fastq")

    @property
    def collapsed_to_hg_rep_sam(self):
        """sorted alignments mapping isoforms to SIRV"""
        return op.join(self.collapse_to_hg_dir, "touse.rep.sam")

    @property
    def collapsed_to_hg_gff(self):
        """gff file of rep isoforms collapsed to hg"""
        return op.join(self.collapse_to_hg_dir, "touse.gff")

    @property
    def collapsed_to_sirv_gff(self):
        """gff file of rep isoforms collapsed to sirv"""
        return op.join(self.collapse_to_sirv_dir, "touse.gff")

    @property
    def collapsed_to_hg_abundance(self):
        """abundance file of rep isoforms collapsed to hg"""
        return op.join(self.collapse_to_hg_dir, "touse.abundance.txt")

    @property
    def collapsed_to_sirv_abundance(self):
        """abundance file of rep isoforms collapsed to sirv"""
        return op.join(self.collapse_to_sirv_dir, "touse.abundance.txt")

    @property
    def collapsed_to_hg_group(self):
        """group file of rep isoforms collapsed to hg"""
        return op.join(self.collapse_to_hg_dir, "touse.group.txt")

    @property
    def collapsed_to_sirv_group(self):
        """abundance file of rep isoforms collapsed to sirv"""
        return op.join(self.collapse_to_sirv_dir, "touse.group.txt")

    @property
    def matchAnnot_out(self):
        """MatchAnot output"""
        return self.collapsed_to_hg_rep_sam + ".matchAnnot.txt"

    @property
    def collapsed_to_sirv_rep_sam(self):
        """sorted alignments mapping isoforms to SIRV"""
        return op.join(self.collapse_to_sirv_dir, "touse.rep.sam")

    @property
    def collapsed_to_hg_rep_readlength_csv(self):
        """readlength csv of representative isoforms collapsed to hg"""
        return op.join(self.csv_dir, "collapsed_to_hg_rep_readlength.csv")

    @property
    def collapsed_to_sirv_rep_readlength_csv(self):
        """readlength csv of representative isoforms collapsed to hg"""
        return op.join(self.csv_dir, "collapsed_to_hg_rep_readlength.csv")

    def __str__(self):
        return "IsoSeq validation files under: %s" % self.root_dir


class ValidationRunner(ValidationFiles):
    """Validating from a smrtlink job dir."""

    def __init__(self, root_dir, smrtlink_job_dir):
        super(ValidationRunner, self).__init__(root_dir=root_dir)
        self.smrtlink_job_dir = op.join(smrtlink_job_dir)
        self.subreads_xml = SMRTLinkIsoSeqFiles(self.smrtlink_job_dir).subreads_xml

    @property
    def all_files(self):
        """Return a list of all files and dirs as [(file_description, file_path)]"""
        return self.common_files + self.human_files + self.sirv_files

    @property
    def common_files(self):
        """all common files, not validation files."""
        return [
            ("root_dir", self.root_dir),
            ("smrtlink_dir", self.smrtlink_job_dir),
            ("validation_report_csv", self.validation_report_csv),

            ("polymerase_readlength_csv", self.polymerase_readlength_csv),
            ("ccs_readlength_csv", self.ccs_readlength_csv),
            ("flnc_readlength_csv", self.flnc_readlength_csv),
            ("consensus_isoforms_readlength_csv", self.consensus_isoforms_readlength_csv),
            ("hq_readlength_csv", self.hq_readlength_csv),
            ("lq_readlength_csv", self.lq_readlength_csv),

            ("isoseq_flnc_fasta", self.isoseq_flnc_fa),
            ("isoseq_nfl_fasta", self.isoseq_nfl_fa),
            ("consensus_isoforms_fasta", self.consensus_isoforms_fa),
            ("hq_isoforms_fasta", self.hq_isoforms_fa)
        ]

    @property
    def human_files(self):
        """human related files."""
        return [
            ("gencode_gtf", self.gencode_gtf),
            ("collapsed_to_hg_rep_fastq", self.collapsed_to_hg_rep_fq),
            ("collapsed_to_hg_rep_sam", self.collapsed_to_hg_rep_sam),
            ("collapsed_to_hg_rep_readlength_csv", self.collapsed_to_hg_rep_readlength_csv),
            ("matchAnnot_out", self.matchAnnot_out)
        ]

    @property
    def sirv_files(self):
        """sirv related files."""
        return [
            ("collapsed_to_sirv_rep_fastq", self.collapsed_to_sirv_rep_fq),
            ("collapsed_to_sirv_rep_sam", self.collapsed_to_sirv_rep_sam),
            ("collapsed_to_sirv_rep_readlength_csv", self.collapsed_to_sirv_rep_readlength_csv),
            ("chained_ids_txt", self.chained_ids_txt)
        ]

    def ln_gencode_gtf(self, gencode_gtf):
        """Make link of gencode_gtf"""
        log.info("Making soft link of gencode annotation gtf")
        ln(gencode_gtf, self.gencode_gtf)

    def ln_sirv_truth_dir(self, sirv_truth_dir):
        """Make a link of sirv ground truth dir"""
        log.info("Making soft link of sirv ground truth dir")
        ln(sirv_truth_dir, self.sirv_truth_dir)

    def make_all_files_from_SMRTLink_job(self):
        """Make all data files from a smrtlink job, including
        * consolidate flnc and nfl xml files to fasta
        * collect ccs, classify, cluster reports from sl job and make validation_report_csv
        * make read length csv files
        """
        mkdir(self.root_dir)
        mkdir(self.fasta_dir)
        mkdir(self.csv_dir)
        mkdir(self.chain_sample_dir)

        smrtlink_job_dir = self.smrtlink_job_dir
        self.make_reports_from_SMRTLink_job(smrtlink_job_dir)
        self.consolidate_xml_from_SMRTLink_job(smrtlink_job_dir)
        #symlink smrtlink_job_dir and files to validation dir
        self.ln_files_from_SMRTLink_job(smrtlink_job_dir)
        self.make_readlength_csvs()

    def make_readlength_csvs(self):
        """Make all read length csv files."""
        log.info("make all readlength csv files.")
        z = [
            (self.ccs_fa_gz, self.ccs_readlength_csv),
            (self.isoseq_flnc_fa, self.flnc_readlength_csv),
            (self.isoseq_nfl_fa, self.nfl_readlength_csv),
            (self.hq_isoforms_fa, self.hq_readlength_csv),
            (self.lq_isoforms_fa, self.lq_readlength_csv),
            (self.consensus_isoforms_fa, self.consensus_isoforms_readlength_csv)
        ]
        for fasta_fn, csv_fn in z:
            make_readlength_csv(fasta_fn=fasta_fn, csv_fn=csv_fn)
        self.make_readlength_csv_for_polymerase()

    def make_readlength_csv_for_polymerase(self):
        """Make approximated read length csv for polymerase reads"""
        make_polymerase_readlength_csv(subreads_xml=self.subreads_xml,
                                       csv_fn=self.polymerase_readlength_csv)

    def make_readlength_csv_for_sirv_isoforms(self):
        """Make read length csv for representative isoforms collapsing to sirv"""
        make_readlength_csv(fasta_fn=self.collapsed_to_sirv_rep_fq,
                            csv_fn=self.collapsed_to_sirv_rep_readlength_csv)

    def make_readlength_csv_for_hg_isoforms(self):
        """Make read length csv for representative isoforms collapsing to hg"""
        make_readlength_csv(fasta_fn=self.collapsed_to_hg_rep_fq,
                            csv_fn=self.collapsed_to_hg_rep_readlength_csv)

    def consolidate_xml_from_SMRTLink_job(self, smrtlink_job_dir):
        """Consolidate xml files from SMRTLink job dir"""
        log.info("consolidate xml files from smrtlink job")
        sl_job = SMRTLinkIsoSeqFiles(smrtlink_job_dir)

        # Consolidate isoseq_flnc.fasta
        consolidate_xml(src=sl_job.isoseq_flnc_ds, dst=self.isoseq_flnc_fa)

        # Consolidate isoseq_nfl.fasta
        consolidate_xml(src=sl_job.isoseq_nfl_ds, dst=self.isoseq_flnc_fa)

    def ln_files_from_SMRTLink_job(self, smrtlink_job_dir):
        """Make soft links to existing isoseq output fasta|fastq files."""
        log.info("make soft links from smrtlink job")
        sl_job = SMRTLinkIsoSeqFiles(smrtlink_job_dir)
        # Make a link of smrtlink dir
        ln(smrtlink_job_dir, op.join(self.root_dir, op.basename(sl_job.root_dir)))

        # Make a link of consensus isoforms fa, hq|lq isoforms fa|fq
        ln(src=sl_job.consensus_isoforms_fa, dst=self.consensus_isoforms_fa)
        ln(src=sl_job.hq_isoforms_fa, dst=self.hq_isoforms_fa)
        ln(src=sl_job.hq_isoforms_fq, dst=self.hq_isoforms_fq)
        ln(src=sl_job.lq_isoforms_fa, dst=self.lq_isoforms_fa)
        ln(src=sl_job.lq_isoforms_fq, dst=self.lq_isoforms_fq)

        ln(src=sl_job.ccs_fa_gz, dst=self.ccs_fa_gz)

    def make_reports_from_SMRTLink_job(self, smrtlink_job_dir):
        """Get reports from a SMRTLink job, including ccs report,
        classify report, cluster report, and write to self.validation_report_csv"""
        log.info("make reports from smrtlink job")
        sl_job = SMRTLinkIsoSeqFiles(smrtlink_job_dir)
        reports_fn = [sl_job.ccs_report_json,
                      sl_job.classify_report_json,
                      sl_job.cluster_report_json]
        d = dict()
        for report_fn in reports_fn:
            print report_fn
            for key, val in json_to_attr_dict(report_fn).iteritems():
                d[key] = val

        #write d to validation_report_csv
        with open(self.validation_report_csv, 'w') as writer:
            writer.write("name\tvalue\n")
            for key in sorted(d.keys()):#d.iteritems():
                writer.write("%s\t%s\n" % (key, d[key]))

    def make_readme_txt(self, args, human_only=False, sirv_only=False):
        """Write args and data files to README file."""
        with open(self.readme_txt, 'w') as writer:
            log.info("args=%s\n", args)
            writer.write("# Created by pbtranscript.testkit.ValidationRunner.make_readme_txt()\n")
            writer.write("args=%s\n\n" % args)

            files = self.all_files # by default show all files
            if human_only:
                files = self.common_files + self.human_files
            elif sirv_only:
                files = self.common_files + self.sirv_files

            for desc, fn in files:
                writer.write("%s=%s\n" % (desc, fn))
