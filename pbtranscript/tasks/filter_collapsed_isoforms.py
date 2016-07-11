#!/usr/bin/env python

"""
Filter out collapsed isoforms whose supportive FL reads
is less than min_count and filter out collapsed isoforms
which are a subset of another isoform.

Note that input collapsed isoforms FASTQ file and output
filtered collapsed isoforms FASTQ file must be specified.
Associated files such as abundance.txt, group.txt, gff.txt
will be inferred from FASTQ file.
"""

import sys
import logging

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptOptions import get_base_contract_parser
from pbtranscript.Utils import rmpath, mv, realpath
from pbtranscript.filtering import filter_by_count, filter_out_subsets

import pbtranscript.tasks.collapse_mapped_isoforms as ci


log = logging.getLogger(__name__)


class Constants(object):
    """Constants used in tool contract."""
    TOOL_ID = "pbtranscript.tasks.filter_collapsed_isoforms"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__

    MIN_COUNT_ID = "pbtranscript.task_options.min_fl_count"
    MIN_COUNT_DEFAULT = 2
    MIN_COUNT_DESC = "Minimum FL count to not filter a collapsed isoform"

    FILTER_OUT_SUBSETS_DEFAULT = True
    FILTER_OUT_SUBSETS_DESC = "Filter out collapsed isoforms which are a subset of another isoform"


def add_filter_collapsed_isoforms_io_arguments(arg_parser):
    """Add io arguments to parser."""
    helpstr = "Input collapsed isoforms in FASTQ"
    arg_parser.add_argument("in_rep_fastq", type=str, help=helpstr)

    helpstr = "Output prefix"
    arg_parser.add_argument("out_rep_fastq", type=str, help=helpstr)
    return arg_parser


def add_filter_collapsed_isoforms_arguments(arg_parser):
    """Add option arguments to parser"""
    f_group = arg_parser.add_argument_group("filter collapsed isoforms arguments")
    f_group.add_argument("--min_count", type=int, default=Constants.MIN_COUNT_DEFAULT,
                         help=Constants.MIN_COUNT_DESC)
    f_group.add_argument("--no_filter_subsets", dest="filter_out_subsets",
                         default=Constants.FILTER_OUT_SUBSETS_DEFAULT,
                         action="store_false", help=Constants.FILTER_OUT_SUBSETS_DESC)
    return arg_parser


def args_runner(args):
    """Run given input args, e.g.,
    filter_collapsed_isoforms.py in_rep_fastq out_rep_fastq --min_count 2
    filter_collapsed_isoforms.py in_rep_fastq out_rep_fastq --min_count 2 --no_filter_subsets
    """
    in_fq, out_fq = args.in_rep_fastq, args.out_rep_fastq

    def _get_prefix_of_rep_fq(fn):
        """Return prefix of *.rep.fq"""
        if fn.endswith(".rep.fastq") or fn.endswith(".rep.fq"):
            return '.'.join(fn.split(".")[0:-2])
        elif fn.endswith(".fastq") or fn.endswith(".fq"):
            return '.'.join(fn.split(".")[0:-1])
        raise ValueError("Invalid collapsed isoforms .rep.fastq file %s" % fn)

    input_prefix = _get_prefix_of_rep_fq(in_fq)
    output_prefix = _get_prefix_of_rep_fq(out_fq)

    # infer group.txt, abundance.txt and gff
    in_group_filename = input_prefix + ".group.txt"
    in_abundance_filename = input_prefix + ".abundance.txt"
    in_gff_filename = input_prefix + ".gff"

    tmp_out_abundance_filename = output_prefix + ".has_subsets.abundance.txt"
    tmp_out_gff_filename = output_prefix + ".has_subsets.gff"
    tmp_out_fq = output_prefix + ".has_subsets.rep.fastq"

    out_abundance_filename = output_prefix + ".abundance.txt"
    out_gff_filename = output_prefix + ".gff"

    # Filter collapsed isoforms by min FL count.
    logging.info("Filtering collapsed isoforms by count %s", args.min_count)
    filter_by_count(in_group_filename=in_group_filename,
                    in_abundance_filename=in_abundance_filename,
                    in_gff_filename=in_gff_filename, in_rep_filename=in_fq,
                    out_abundance_filename=tmp_out_abundance_filename,
                    out_gff_filename=tmp_out_gff_filename, out_rep_filename=tmp_out_fq,
                    min_count=args.min_count)

    # Remove collapsed isoforms which are a subset of another isoform
    logging.info("Filtering out subsets collapsed isoforms = %s", args.filter_out_subsets)
    if args.filter_out_subsets is True:
        filter_out_subsets(in_abundance_filename=tmp_out_abundance_filename,
                           in_gff_filename=tmp_out_gff_filename,
                           in_rep_filename=tmp_out_fq,
                           out_abundance_filename=out_abundance_filename,
                           out_gff_filename=out_gff_filename,
                           out_rep_filename=out_fq,
                           max_fuzzy_junction=args.max_fuzzy_junction)
        rmpath(tmp_out_abundance_filename)
        rmpath(tmp_out_gff_filename)
        rmpath(tmp_out_fq)
    else:
        mv(tmp_out_abundance_filename, out_abundance_filename)
        mv(tmp_out_gff_filename, out_gff_filename)
        mv(tmp_out_fq, out_fq)

    logging.info("Filtered collapsed isoforms sequences written to %s", realpath(out_fq))
    logging.info("Filtered collapsed isoforms abundance written to %s", realpath(out_abundance_filename))
    logging.info("Filtered collapsed isoforms gff written to %s", realpath(out_gff_filename))
    return 0


def resolved_tool_contract_runner(rtc):
    """Run given a resolved tool contract"""
    raise NotImplementedError() # Merged to post_mapping_to_genome


def add_filter_collapsed_isoforms_tcp_options(tcp):
    """Add tcp options"""
    tcp.add_int(option_id=Constants.MIN_COUNT_ID, option_str="min FL count",
                name="min FL count", default=Constants.MIN_COUNT_DEFAULT,
                description=Constants.MIN_COUNT_DESC)
    tcp.add_int(option_id=ci.Constants.MAX_FUZZY_JUNCTION_ID,
                option_str="max_fuzzy_junction", name="max fuzzy junction",
                default=ci.Constants.MAX_FUZZY_JUNCTION_DEFAULT,
                description=ci.Constants.MAX_FUZZY_JUNCTION_DESC)
    return tcp


def get_contract_parser():
    """Get tool contract parser.
    Input:
        idx 0 - group txt file
        idx 1 - abundance txt file
        idx 2 - gff file
        idx 3 - rep fastq file
    Output:
        idx 0 - abundance file
        idx 1 - gff file
        idx 2 - rep fastq file
    """
    p = get_base_contract_parser(Constants, default_level="DEBUG")

    # argument parser
    add_filter_collapsed_isoforms_io_arguments(p.arg_parser.parser)
    add_filter_collapsed_isoforms_arguments(p.arg_parser.parser)
    # DONT move --max_fuzzy_junction to add_filter_collapsed_isoforms_arguments
    p.arg_parser.parser.add_argument("--max_fuzzy_junction", type=int,
                                     default=ci.Constants.MAX_FUZZY_JUNCTION_DEFAULT,
                                     help=ci.Constants.MAX_FUZZY_JUNCTION_DESC)

    # tool contract parser
    tcp = p.tool_contract_parser
    tcp.add_input_file_type(FileTypes.TXT, "group_txt", "TXT In",
                            "Group file associating isoforms with reads") # input 0
    tcp.add_input_file_type(FileTypes.TXT, "abundance_txt", "TXT In",
                            "Abundance file") # input 1
    tcp.add_input_file_type(FileTypes.GFF, "collapsed_isoforms_gff", "GFF In",
                            "GFF file") # input 2
    tcp.add_input_file_type(FileTypes.FASTQ, "rep_fastq", "FASTQ In",
                            "Representative sequences of collapsed isoforms in FATSQ") # input 3

    tcp.add_output_file_type(FileTypes.TXT, "abundance_txt",
                             name="TXT file", description="Abundance file",
                             default_name="output_mapped_filtered_abundance") # output 0
    tcp.add_output_file_type(FileTypes.GFF, "collapsed_isoforms_gff",
                             name="GFF file", description="Collapsed isoforms in GFF file",
                             default_name="output_mapped_filtered") # output 1
    tcp.add_output_file_type(FileTypes.FASTQ, "rep_fastq", name="FASTQ file",
                             description="Representative sequences of collapsed isoforms in FASTQ",
                             default_name="output_mapped_filtered") # output 2
    add_filter_collapsed_isoforms_tcp_options(tcp)
    return p


def main(args=sys.argv[1:]):
    """Main"""
    return pbparser_runner(argv=args,
                           parser=get_contract_parser(),
                           args_runner_func=args_runner,
                           contract_runner_func=resolved_tool_contract_runner,
                           alog=log,
                           setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main())
