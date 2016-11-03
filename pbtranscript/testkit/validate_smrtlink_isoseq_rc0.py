#!/usr/bin/env python

"""
Validate a smrtlink RC0 isoseq job with reference and gencode.
"""


import argparse
import sys
import logging
import os.path as op

from pbtranscript.Utils import execute, realpath, ln, mkdir
from pbtranscript.tasks.post_mapping_to_genome import add_post_mapping_to_genome_arguments, \
                post_mapping_to_genome_runner
from pbtranscript.collapsing.CollapsingUtils import map_isoforms_and_sort
from pbtranscript.tasks.map_isoforms_to_genome import add_gmap_arguments


__author__ = 'etseng@pacb.com, yli@pacb.com'

FORMATTER = op.basename(__file__) + ':%(levelname)s:'+'%(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMATTER)
log = logging.getLogger(__name__)


GMAP_DB = "/pbi/dept/secondary/siv/testdata/isoseq/gmap_db/"
GMAP_NAME = "hg38"
GMAP_NPROC = 12
GENCODE_GTF = "/pbi/dept/secondary/siv/testdata/isoseq/gencode/gencode.v25.annotation.gtf"


def smrtlink_file(smrtlink_job_dir, basename):
    """
    Return {smrtlink_job_dir}/tasks/pbtranscript.tasks.combine_cluster_bins-0/basename
    """
    return op.join(smrtlink_job_dir, 'tasks', 'pbtranscript.tasks.combine_cluster_bins-0', basename)


def link_files(smrtlink_job_dir, out_dir, more_files):
    """
    Make soft link of some smrtlink isoseq job output files and more_files in {out_dir}.
    """
    log.info("Making soft link of files")
    hq_fq = smrtlink_file(smrtlink_job_dir=smrtlink_job_dir, basename="hq_isoforms.fastq")
    cluster_report = smrtlink_file(smrtlink_job_dir=smrtlink_job_dir, basename="cluster_report.csv")
    hq_lq_prefix_pickle = smrtlink_file(smrtlink_job_dir=smrtlink_job_dir,
                                        basename="hq_lq_prefix_dict.pickle")

    assert isinstance(more_files, list)
    fs = more_files + [hq_fq, cluster_report, hq_lq_prefix_pickle]
    for f in fs:
        dst = op.join(out_dir, op.basename(f))
        log.debug("%s --> %s", f, dst)
        ln(f, dst)

    return hq_fq, hq_lq_prefix_pickle


def collapse_to_reference(hq_fq, hq_lq_prefix_pickle, out_dir, args, out_prefix="touse"):
    """
    hq_fq --- HQ isoforms in fastq format.
    out_dir --- where to output results.
    min_count --- minimum # of supportive FLNC reads to call an isoform HQ

    First map HQ isoforms against reference ($gmap_db_dir/$gmap_db_name),
    next collapse HQ isoforms to representative isoforms based on mapping,
    and finally map representative isoforms to reference and return sorted
    SAM output.
    """
    gmap_db, gmap_name = args.gmap_db, args.gmap_name
    # Map HQ isoforms to GMAP reference genome
    log.info("Mapping HQ isoforms to reference.")
    log.debug("HQ isoforms: %s", hq_fq)
    log.debug("reference: %s/%s", gmap_db, gmap_name)

    hq_sam = op.join(out_dir, "%s.sorted.sam" % op.basename(hq_fq))
    map_isoforms_and_sort(input_filename=hq_fq, sam_filename=hq_sam,
                          gmap_db_dir=gmap_db, gmap_db_name=gmap_name, gmap_nproc=GMAP_NPROC)

    log.info("Collapsing and filtering HQ isoforms to create representative isoforms.")
    # Post mapping to genome analysis, including
    #   * collapse polished HQ isoform clusters into groups
    #   * count abundance of collapsed isoform groups
    #   * filter collapsed isoforms based on abundance info
    rep_fq = op.join(out_dir, "%s.rep.fastq" % out_prefix)
    post_mapping_to_genome_runner(in_isoforms=hq_fq, in_sam=hq_sam,
                                  in_pickle=hq_lq_prefix_pickle,
                                  out_isoforms=rep_fq, out_gff=None, out_abundance=None,
                                  out_group=None, out_read_stat=None,
                                  min_aln_coverage=args.min_aln_coverage,
                                  min_aln_identity=args.min_aln_identity,
                                  min_flnc_coverage=args.min_flnc_coverage,
                                  max_fuzzy_junction=args.max_fuzzy_junction,
                                  allow_extra_5exon=args.allow_extra_5exon,
                                  min_count=args.min_count)

    rep_sam = op.join(out_dir, "%s.sorted.sam" % rep_fq)
    # Map representitive isoforms to reference
    map_isoforms_and_sort(input_filename=rep_fq, sam_filename=rep_sam,
                          gmap_db_dir=gmap_db, gmap_db_name=gmap_name, gmap_nproc=GMAP_NPROC)
    return rep_sam


def validate_with_Gencode(sorted_rep_sam, gencode_gtf, match_out):
    """
    Input:
      sorted_rep_sam -- sorted SAM output mapping (collapsed) representitve isoforms to reference
      eval_dir -- evaluation directory
    Run matchAnnot to compare sorted_rep_sam with gencode v25 and output to eval_dir
    """
    log.info("Writing matchAnnot output to %s", match_out)
    cmd = "matchAnnot.py --gtf={0} {1} > {2}".format(gencode_gtf, sorted_rep_sam, match_out)
    execute(cmd)


def check_matchAnnot_out(match_out, min_percentage):
    """
    Check from matchAnnot.txt % of isoforms with scores are 4 or 5
    """
    lines = [l for l in open(match_out, 'r') if l.startswith("summary:")]
    total_n = 0
    ns = [0, 0, 0, 0, 0, 0] # ns[i] number of isoforms scores == i
    for l in lines:
        if "isoforms read" in l:
            total_n = [int(s) for s in l.split() if s.isdigit()][0]

        for score in [5, 4, 3, 2, 1, 0]:
            if "isoforms scored %s" % score in l:
                ns[score] = [int(s.replace(',', '')) for s in l.split()
                             if s.replace(',', '').isdigit()][0]
    if sum(ns) != total_n:
        raise ValueError("MatchAnnot %s (%s isoforms read != %s isoforms classified)"
                         % (match_out, total_n, sum(ns)))

    if total_n * min_percentage / 100.0 > ns[5] + ns[4]:
        # isoforms with score are 4 or 5 is less than min_percentage
        log.info("There are %s out of %s isoforms with scores equal 4 or 5 < %s percent.",
                 ns[5]+ns[4], total_n, min_percentage)
        return False
    else:
        log.info("There are %s out of %s isoforms with scores equal 4 or 5 >= %s percent.",
                 ns[5]+ns[4], total_n, min_percentage)
        return True


def make_sane(args):
    """Make sane of input output"""
    args.smrtlink_job_dir = realpath(args.smrtlink_job_dir)
    args.out_dir = realpath(args.out_dir)

    if args.gmap_db is None:
        args.gmap_db = realpath(GMAP_DB)
        log.warning("Reset GMAP DB to %s", args.gmap_db)

    if args.gmap_name is None:
        args.gmap_name = GMAP_NAME
        log.warning("Reset GMAP NAME to %s", args.gmap_name)

    if not op.exists(args.smrtlink_job_dir):
        raise IOError("SMRTLink job directory %s does not exist" % args.smrtlink_job_dir)

    if not op.exists(op.join(args.gmap_db, args.gmap_name)):
        raise IOError("GMAP reference %s/%s does not exist." % (args.gmap_db, args.gmap_name))

    if not op.exists(args.gencode_gtf):
        raise IOError("Gencode gtf file %s does not exist." % args.gencode_gtf)

    log.info("Making out_dir %s", args.out_dir)
    mkdir(args.out_dir)
    return args


def run(args):
    """
    Collapse HQ isoforms from SMRTLink Iso-Seq (w/wo genome) job to hg38.
    """
    args = make_sane(args)

    hq_fq, hq_lq_prefix_pickle = link_files(smrtlink_job_dir=args.smrtlink_job_dir,
                                            out_dir=args.out_dir, more_files=[args.gencode_gtf])

    # Collapse HQ isoforms (hq_fq) to representive isoforms, then map
    # representative isoforms to gmap reference, and sort output SAM (sorted_rep_sam).
    sorted_rep_sam = collapse_to_reference(hq_fq=hq_fq, hq_lq_prefix_pickle=hq_lq_prefix_pickle,
                                           out_dir=args.out_dir, args=args)

    # Run matchAnnot.py to compare sorted_rep_sam against gencode gtf.
    match_out = sorted_rep_sam + ".matchAnnot.txt"
    validate_with_Gencode(sorted_rep_sam=sorted_rep_sam, gencode_gtf=args.gencode_gtf,
                          match_out=match_out)

    # Check from matchAnnot.txt % of isoforms with scores 4 or 5
    res = check_matchAnnot_out(match_out, args.min_percentage)

    if res:
        log.info("Test PASSED.")
    else:
        log.info("Test FAILED.")
        raise ValueError("Test FAILED.")


def get_parser():
    """Get argument parser."""
    helpstr = "Validate SMRTLink Iso-Seq (w/wo Genome) output against GenCode"
    parser = argparse.ArgumentParser(helpstr)

    helpstr = "Smrtlink Iso-Seq (w/wo Genome) job directory"
    parser.add_argument("smrtlink_job_dir", help=helpstr)

    helpstr = "Validation output directory"
    parser.add_argument("out_dir", help=helpstr)

    helpstr = "min percentage of isoforms whose matchAnnot scores are 4 or 5 (default: 80)"
    parser.add_argument("--min_percentage", type=int, default=80, help=helpstr)

    parser = add_post_mapping_to_genome_arguments(parser)
    parser = add_gmap_arguments(parser)

    helpstr = "Gencode gtf file"
    parser.add_argument("--gencode_gtf", type=str, default=GENCODE_GTF, help=helpstr)

    return parser


if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))