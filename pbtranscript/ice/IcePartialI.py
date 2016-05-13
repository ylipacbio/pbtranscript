#!/usr/bin/env python
###############################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

"""
Overview:
    pbtranscript cluster contains two main components:
    * (1) ICE (iterative clustering and error correction) to predict
      unpolished consensus isoforms.
    * (2) Polish, to use nfl reads and quiver to polish those predicted
      unpolished isoforms. Polish contains three steps:
      + (2.1) IceAllPartials (ice_partial.py all)
              Align and assign nfl reads to unpolished isoforms, and
              save results to a pickle file.
      + (2.2) IceQuiver (ice_quiver.py all)
              Call quiver to polish each isoform based on alignments
              created by mapping its associated fl and nfl reads to
              this isoform.
      + (2.3) IceQuiverPostprocess (ice_quiver.py postprocess)
              Collect and post process quiver results, and classify
              HQ/LQ isoforms.

    In order to handle subtasks by SMRTPipe instead of pbtranscript
    itself, we will refactor the polish phase including
    (2.1) (2.2) and (2.3).

    (2.1) IceAllPartials (ice_partial.py all) will be refactored to
      + (2.1.1) ice_partial.py split
                Split nfl reads into N chunks (N<=100).
      + (2.1.2) ice_partial.py i
                For each chunk i, align and assign its reads to unpolished
                isoforms and create a pickle file.
      + (2.1.3) ice_partial.py merge
                Merge pickles for all splitted chunks together to a
                big pickle.

    *** Here we are focusing on (2.1.2) ice_partial.py i ***

Description:
    (2.1.2) ice_partial.py i

    Assumption:
     * Phase (1) ICE is done and unpolished isoforms created.
     * Step (2.1.1) is done, and nfl reads have been splitted into N
       chunks.
     * We are processing the i-th chunk of nfl reads.

    Process:
        Given root_dir, and i, aligns reads in the i-th chunk of nfl reads
        to unpolished consensus isoforms and creates an output pickle file.

        The i-th chunk is:
            root_dir/output/map_noFL/input.split_{0:03d}.fasta
        The unpolished isoforms are saved in:
            root_dir/output/final.consensus.fasta
        If available, the suffix array for unpolished isoforms is at:
            root_dir/output/final.consensus.fasta.sa
        The output pickle file for the i-th chunk is:
            *.partial_uc.pickle
        The done file to indicate that this job is successfully completed:
            *.fasta.partial_uc.pickle.DONE
        The number of chunks, N, is not nececssary for this script.

        Log this commnd to a script file for debugging:
            *.partial_uc.sh

    Input:
        Positional:
            root_dir, an output directory for running pbtranscript cluster.
            i, the i-th chunk, 0 <= i < N.
        Optional:
            ccs_fofn, A FOFN of ccs.h5 files (e.g., reads_of_insert.fofn),
            which contain quality values of reads of insert.
            If not given, assume there is no QV information available.
        blasr_nproc: Number of blasr nproc.

    Output:
        Align reads in the i-th chunk to unpolished consensus isoforms
        using BLASR, assign reads to isoforms and then save the following
        to a pickle:
            {isoform_id:[read_ids], nohit: set(no_hit_read_ids)}
        touch a done file when this job is finished.

    Hierarchy:
        pbtranscript = iceiterative

        pbtranscript --quiver = iceiterative + \
                                ice_polish.py

        ice_polish.py =  ice_make_fasta_fofn.py + \
                         ice_partial.py all + \
                         ice_quiver.py all

        ice_partial.py all = ice_partial.py split + \
                             ice_partial.py i + \
                             ice_partial.py merge

        (ice_partial.py one --> only apply ice_partial on a given input fasta)

        ice_quiver.py all = ice_quiver.py i + \
                            ice_quiver.py merge + \
                            ice_quiver.py postprocess

    Example:
        ice_partial.py i root_dir {i} --ccs_fofn=ccs_fofn \
                         --blasr_nproc=blasr_nproc

"""

import logging
from pbtranscript.__init__ import get_version
from pbtranscript.PBTranscriptOptions import add_fofn_arguments, \
    add_cluster_root_dir_as_positional_argument, add_tmp_dir_argument
from pbtranscript.Utils import nfs_exists
from pbtranscript.ice.IceFiles import IceFiles
from pbtranscript.ice.IcePartial import build_uc_from_partial_daligner
from pbtranscript.ice.__init__ import ICE_PARTIAL_PY


def add_ice_partial_i_arguments(parser):
    """Add IcePartialI arguments."""
    parser = add_cluster_root_dir_as_positional_argument(parser)

    helpstr = "To process the i-th chunk of non-full-length reads."
    parser.add_argument("i", nargs="+", help=helpstr, type=int)

    parser.add_argument("--blasr_nproc",
                        type=int,
                        dest="blasr_nproc",
                        action="store",
                        default=12,
                        help="Number of cores for each BLASR job.")

    parser = add_fofn_arguments(parser, ccs_fofn=True)
    parser = add_tmp_dir_argument(parser)
    return parser


class IcePartialI(object):

    """IcePartialI object."""

    # Define description of IcePartialI
    desc = "Process the i-th chunk of non-full-length reads, align " + \
           "and assign these reads to unpolished consensus isoforms " + \
           "and create a pickle file."

    prog = "%s i " % ICE_PARTIAL_PY  # used by cmd_str

    def __init__(self, root_dir, i, ccs_fofn, blasr_nproc, tmp_dir):
        """
        root_dir --- root directory for saving intermediate files running
        pbtranscript cluster.
        i --- a list of nfl 1 or more read chunks to assign to isoforms
        ccs_fofn --- FOFN of ccs.h5
        blasr_nproc --- blasr nproc
        tmp_dir - where to save intermediate files such as dazz files.
                  if None, writer dazz files to the same directory as
                  input query/target.
        """
        self.root_dir = root_dir
        self.i = i
        assert(isinstance(self.i, list))
        self.ccs_fofn = ccs_fofn
        self.blasr_nproc = int(blasr_nproc)
        self.tmp_dir = tmp_dir

    def getVersion(self):
        """Return version string."""
        return get_version()

    def cmd_str(self):
        """Return a cmd string of this object."""
        return self._cmd_str(root_dir=self.root_dir, i=self.i,
                             ccs_fofn=self.ccs_fofn,
                             blasr_nproc=self.blasr_nproc,
                             tmp_dir=self.tmp_dir)

    def _cmd_str(self, root_dir, i, ccs_fofn, blasr_nproc, tmp_dir):
        """Return a cmd string given parameters."""
        cmd = self.prog + \
            "{d} ".format(d=root_dir) + \
            "{i} ".format(i=" ".join([str(x) for x in i]))
        if ccs_fofn is not None:
            cmd += "--ccs_fofn={ccs_fofn} ".format(ccs_fofn=ccs_fofn)
        if blasr_nproc is not None:
            cmd += "--blasr_nproc={n} ".format(n=blasr_nproc)
        if tmp_dir is not None:
            cmd += "--tmp_dir={t} ".format(t=tmp_dir)
        return cmd

    def _validate_inputs(self, root_dir, i, ccs_fofn, blasr_nproc, tmp_dir):
        """
        Check inputs, write $ICE_PARTIAL_PY i command to script_file
        and return (input_fasta, ref_fasta, out_pickle, done_file)
        for the i-th chunk of nfl reads.
        """
        icef = IceFiles(prog_name="ice_partial_{i}".format(i=i),
                        root_dir=root_dir, no_log_f=False)

        # root_dir/output/final.consensus.fasta
        ref_fasta = icef.final_consensus_fa
        ref_dazz = icef.final_dazz_db

        # root_dir/output/map_noFL/input.split_{0:03d}.fasta
        input_fasta = icef.nfl_fa_i(i)

        # $input_fasta.partial_uc.pickle
        out_pickle = icef.nfl_pickle_i(i)

        # $input_fasta.partial_uc.pickle.DONE
        done_file = icef.nfl_done_i(i)

        # $input_fasta.partial_uc.sh
        script_file = icef.nfl_script_i(i)

        # Check if inputs exist.
        errMsg = ""
        if not nfs_exists(input_fasta):
            errMsg = ("The {i}-th splitted non-full-length reads ".format(i=i) +
                      "fasta file {f} does not exist. ".format(f=input_fasta) +
                      "Please run $ICE_PARTIAL_PY split first.")
        elif not nfs_exists(ref_fasta):
            errMsg = ("The unpolished consensus isoforms fasta file " +
                      "{f} does not exist. ".format(f=ref_fasta) +
                      "Please make sure ICE is successfully done.")
        elif not nfs_exists(ref_dazz):
            errMsg = ("The dazz db " +
                      "{f} does not exist. ".format(f=ref_dazz) +
                      "Please make sure it is already built.")
        if len(errMsg) != 0:
            raise ValueError(errMsg)

        # Save cmd to script_file.
        cmd = self._cmd_str(root_dir=root_dir, i=[i],
                            ccs_fofn=ccs_fofn,
                            blasr_nproc=blasr_nproc,
                            tmp_dir=tmp_dir)
        with open(script_file, 'w') as writer:
            writer.write(cmd + "\n")

        icef.add_log("Writing CMD to: {script_file}".
                     format(script_file=script_file))
        icef.close_log()

        return (input_fasta, ref_fasta, out_pickle, done_file)

    def run(self):
        """Run IcePartialI"""
        logging.debug("root_dir: {d}.".format(d=self.root_dir))
        logging.debug("i: {f}.".format(f=self.i))
        logging.debug("ccs_fofn: {ccs_fofn}.".format(ccs_fofn=self.ccs_fofn))

        # Validate input files, write equivalent command
        # to script_file.
        for i in self.i:
            input_fasta, ref_fasta, out_pickle, done_file = \
                self._validate_inputs(root_dir=self.root_dir,
                                      i=i, ccs_fofn=self.ccs_fofn,
                                      blasr_nproc=self.blasr_nproc,
                                      tmp_dir=self.tmp_dir)

            logging.info("Calling build_uc_from_partial_daligner.")
            build_uc_from_partial_daligner(input_fasta=input_fasta,
                                           ref_fasta=ref_fasta,
                                           out_pickle=out_pickle,
                                           tmp_dir=self.tmp_dir,
                                           ccs_fofn=self.ccs_fofn,
                                           done_filename=done_file,
                                           cpus=self.blasr_nproc,
                                           no_qv_or_aln_checking=True)
