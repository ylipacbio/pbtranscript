#!/usr/bin/env python
import os.path as op
import sys
import logging
from argparse import ArgumentParser
from cPickle import load, dump
from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbtranscript.ice.ProbModel import ProbFromModel, ProbFromFastq
import pbtranscript.ice.IceIterative as ice
from pbtranscript.ice.IceUtils import ice_fa2fq, realpath
from pbtranscript.PBTranscriptOptions import add_fofn_arguments

FORMATTER = op.basename(__file__) + ':%(levelname)s:'+'%(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMATTER)
log = logging.getLogger(__name__)


__author__ = "etseng@pacificbiosciences.com"


def add_picking_up_ice_arguments(parser):
    """Add arguments for picking up ice."""
    helpstr = "Last successful pickle (e.g., clusterOut/output/input.split_001.fa.pickle)"
    parser.add_argument("pickle_filename", help=helpstr)

    helpstr = "root dir (default: cluster_out/)"
    parser.add_argument("--root_dir", default="cluster_out", help=helpstr)

    # FIXME, use ccs.xml instead
    parser = add_fofn_arguments(parser, ccs_fofn=True)

    helpstr = "Input full-length non-chimeric reads in FASTA or ContigSet format, " + \
              "used for clustering consensus isoforms, e.g., isoseq_flnc.fasta"
    parser.add_argument("--flnc", "--flnc_fa", default="isoseq_flnc.fasta", help=helpstr)

    helpstr = "Comma-separated additional fasta files to add."
    parser.add_argument("--fasta_files_to_add", default=None, help=helpstr)
    return parser


def current_fasta(root_dir):
    """Return filename of root_dir/current.fasta, which is used in
    the `current` ICE iteration."""
    return realpath(op.join(root_dir, "current.fasta"))


def current_fastq(root_dir):
    """Return filename of root_dir/current.fastq, which is used in
    the `current` ICE iteration."""
    return realpath(op.join(root_dir, "current.fastq"))


def ensure_pickle_goodness(pickle_filename, root_dir, fasta_files_to_add=None):
    """
    Old versions of IceIterative.write_pickle is missing some key/values.
    Add if needed.
    Return a good pickle object, and the pickle's filename
    """
    a = load(open(pickle_filename))
    if realpath(a['fasta_filename']) != current_fasta(root_dir):
        raise ValueError("The pickle file %s indicates that " % pickle_filename +
                         "current.fasta is not being used. ICE likely did not " +
                         "finish to a point that could be picked up.")

    a['newids'] = check_n_fix_newids(a)

    if fasta_files_to_add is not None:
        for f in fasta_files_to_add.split(','):
            if not op.exists(f):
                raise IOError("%s is not a valid fasta file to add!" % f)
            if f in a['fasta_filenames_to_add']:
                log.warning("%s is already in to-add list. Ignore.", f)
            a['fasta_filenames_to_add'].append(f)
    if 'root_dir' not in a:
        log.warning("Pickle %s missing some key-values. Fixing it.", pickle_filename)
        a['root_dir'] = root_dir
        a['all_fasta_filename'] = a['all_fasta_fiilename']
        a['qv_prob_threshold'] = 0.03
        fixed_pickle_filename = pickle_filename + ".fixed"
        with open(fixed_pickle_filename, 'w') as f:
            dump(a, f)
        log.info("Fixed pickle written to %s", fixed_pickle_filename)
        return a, fixed_pickle_filename
    else:
        # newid might have been fixed, STILL output pickle writing anyway
        with open(pickle_filename, 'w') as f:
            dump(a, f)
        return a, pickle_filename


def check_n_fix_newids(icec_obj):
    """
    Check and fix newids in icec_obj['newids'], return valid newids
    """
    newids = icec_obj['newids']

    if len(newids) == 0:
        log.info("newids is empty (probably a finalized run). set it.")
        for k, v in icec_obj['d'].iteritems():
            if len(v) != 1:
                newids.add(k)
        log.info("added %s seqs to newids", len(newids))
    return newids


def make_current_fasta(icec_obj, flnc_filename, root_dir):
    """
    current fasta will consists of all ids

    however --- if this was a already finished run and we are adding more input,
        then newids is empty, in this case we set newids = everything that
        has no affiliation or more than one affiliated cluster in d
    """
    with FastaWriter(current_fasta(root_dir)) as f:
        for r in FastaReader(flnc_filename):
            f.writeRecord(r)


def pickup_icec_job(pickle_filename, ccs_fofn, flnc_filename, fasta_files_to_add, root_dir):
    """
    Reconstruct an ICE object from a pickle file and restart this ICE job.
    """
    log.info("Reading ICE pickle %s ....", pickle_filename)
    icec_obj, icec_pickle_filename = ensure_pickle_goodness(
        pickle_filename=pickle_filename, root_dir=root_dir,
        fasta_files_to_add=fasta_files_to_add)

    c_fa = current_fasta(root_dir)
    c_fq = current_fastq(root_dir)
    log.info("Making current.fasta %s for ICE ....", c_fa)
    make_current_fasta(icec_obj=icec_obj, flnc_filename=flnc_filename, root_dir=root_dir)

    log.info("Loading prob QV information....")
    probqv = None
    if ccs_fofn is None:
        logging.info("Loading probability from model (0.01,0.07,0.06)")
        probqv = ProbFromModel(.01, .07, .06)
    else:
        #if use_finer_qv:
        #    probqv = ProbFromQV(input_fofn=ccs_fofn, fasta_filename=input_fasta)
        #    logging.info("Loading prob QVs from %s + %s took %s secs",
        #                 ccs_fofn, input_fasta, time.time()-start_t)
        logging.info("Converting %s to %s", c_fa, c_fq)
        ice_fa2fq(c_fa, ccs_fofn, c_fq)

        logging.info("Loading prob QVs from %s", c_fq)
        probqv = ProbFromFastq(c_fq)

    log.info("Starting ICE from pickle %s....", icec_pickle_filename)
    icec = ice.IceIterative.from_pickle(icec_pickle_filename, probqv)

    # first must RE-RUN gcon to get all the proper refs
    icec.changes = set()
    icec.refs = {}
    icec.ccs_fofn = ccs_fofn
    icec.all_fasta_filename = flnc_filename
    todo = icec.uc.keys()
    log.info("Re-run gcon for proper refs....")
    icec.run_gcon_parallel(todo)

    log.info("Re-calculating cluster prob, just to be safe....")
    icec.calc_cluster_prob(True)

    log.info("Sanity checking uc_refs now....")
    icec.sanity_check_uc_refs()

    log.info("Ensuring prob QV of new ids are consistent....")
    icec.ensure_probQV_newid_consistency()

    log.info("Sanity check done. Resuming ICE job.")
    icec.run()


def main():
    """Main function to pick up ICE"""
    parser = add_picking_up_ice_arguments(parser=ArgumentParser())
    args = parser.parse_args()
    pickup_icec_job(pickle_filename=args.pickle_filename,
                    ccs_fofn=args.ccs_fofn,
                    flnc_filename=args.flnc,
                    fasta_files_to_add=args.fasta_files_to_add,
                    root_dir=args.root_dir)
    return 0

if __name__ == "__main__":
    sys.exit(main())
