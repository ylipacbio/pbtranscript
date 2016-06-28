#!/usr/bin/env python

"""
Filter collapsed isoforms by count (abundance info).
"""

import sys
import logging

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptOptions import get_base_contract_parser
from pbtranscript.filtering import filter_by_count


log = logging.getLogger(__name__)


class Constants(object):
    """Constants used in tool contract."""
    TOOL_ID = "pbtranscript.tasks.filter_by_count"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__

    MIN_COUNT_ID = "pbtranscript.task_options.min_fl_count"
    MIN_COUNT_DEFAULT = 2
    MIN_COUNT_DESC = "Minimum FL count to not filter a collapsed isoform"


def add_filter_by_count_io_arguments(arg_parser):
    """Add io arguments to parser."""
    helpstr = "Input prefix"
    arg_parser.add_argument("input_prefix", type=str, help=helpstr)

    helpstr = "Output prefix"
    arg_parser.add_argument("output_prefix", type=str, help=helpstr)
    return arg_parser


def add_filter_by_count_arguments(arg_parser):
    """Add option arguments to parser"""
    fbc_group = arg_parser.add_argument_group("filter by count arguments")
    helpstr = "Minimum FL count to consider a collapsed isoform good (default 2)."
    fbc_group.add_argument("--min_count", type=int,
                           default=Constants.MIN_COUNT_DEFAULT, help=helpstr)
    return arg_parser


def args_runner(args):
    """Run given input args, e.g.,
    filter_by_count.py input_prefix output_prefix --min_count 2
    """
    input_prefix, output_prefix = args.prefix, args.output_prefix
    in_group_filename = input_prefix + ".group.txt"
    in_abundance_filename = input_prefix + ".abundance.txt"
    in_gff_filename = input_prefix + ".gff"
    in_rep_filename = input_prefix + ".rep.fastq"

    out_abundance_filename = output_prefix + ".abundance.txt"
    out_gff_filename = output_prefix + ".gff"
    out_rep_filename = output_prefix + ".rep.fastq"

    filter_by_count(in_group_filename=in_group_filename,
                    in_abundance_filename=in_abundance_filename,
                    in_gff_filename=in_gff_filename,
                    in_rep_filename=in_rep_filename,
                    out_abundance_filename=out_abundance_filename,
                    out_gff_filename=out_gff_filename,
                    out_rep_filename=out_rep_filename,
                    min_count=Constants.MIN_COUNT_DEFAULT)
    return 0


def resolved_tool_contract_runner(rtc):
    """Run given a resolved tool contract"""
    raise NotImplementedError() # Merged to post_mapping_to_genome
#    filter_by_count(in_group_filename=rtc.task.input_files[0],
#                    in_abundance_filename=rtc.task.input_files[1],
#                    in_gff_filename=rtc.task.input_files[2],
#                    in_rep_filename=rtc.task.input_files[3],
#                    out_abundance_filename=rtc.task.output_files[0],
#                    out_gff_filename=rtc.task.output_files[1],
#                    out_rep_filename=rtc.task.output_files[2],
#                    min_count=Constants.MIN_COUNT_DEFAULT)
#    return 0


def add_filter_by_count_tcp_options(tcp):
    """Add tcp options"""
    tcp.add_int(option_id=Constants.MIN_COUNT_ID, option_str="min FL count",
                name="min FL count", default=Constants.MIN_COUNT_DEFAULT,
                description=Constants.MIN_COUNT_DESC)
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
    add_filter_by_count_io_arguments(p.arg_parser.parser)
    add_filter_by_count_arguments(p.arg_parser.parser)

    # tool contract parser
#    tcp = p.tool_contract_parser
#    tcp.add_input_file_type(FileTypes.TXT, "group_txt", "TXT In",
#                            "Group file associating isoforms with reads") # input 0
#    tcp.add_input_file_type(FileTypes.TXT, "abundance_txt", "TXT In",
#                            "Abundance file") # input 1
#    tcp.add_input_file_type(FileTypes.GFF, "collapsed_isoforms_gff", "GFF In",
#                            "GFF file") # input 2
#    tcp.add_input_file_type(FileTypes.FASTQ, "rep_fastq", "FASTQ In",
#                            "Representative sequences of collapsed isoforms in FATSQ") # input 3
#
#    tcp.add_output_file_type(FileTypes.TXT, "abundance_txt",
#                             name="TXT file", description="Abundance file",
#                             default_name="output_mapped_filtered_abundance") # output 0
#    tcp.add_output_file_type(FileTypes.GFF, "collapsed_isoforms_gff",
#                             name="GFF file", description="Collapsed isoforms in GFF file",
#                             default_name="output_mapped_filtered") # output 1
#    tcp.add_output_file_type(FileTypes.FASTQ, "rep_fastq", name="FASTQ file",
#                             description="Representative sequences of collapsed isoforms in FASTQ",
#                             default_name="output_mapped_filtered") # output 2
#    add_filter_by_count_tcp_options(tcp)
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
