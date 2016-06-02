#!/usr/bin/env python

"""
Map isoforms to reference genomes and sort.
"""
import sys
import logging

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptOptions import get_base_contract_parser
from pbtranscript.collapsing.CollapsingUtils import map_isoforms_and_sort

log = logging.getLogger(__name__)

class Constants(object):
    """Constants used in tool contract."""
    TOOL_ID = "pbtranscript.tasks.map_isoforms"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__

    GMAP_NAME_ID = "pbtranscript.task_options.gmap_name"
    GMAP_NAME_DEFAULT = "SIRV"

    GMAP_DB_ID = "pbtranscript.task_options.gmap_db"
    GMAP_DB_DEFAULT = "/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data/gmap_db"

    GMAP_NPROC_ID = "pbtranscript.task_options.gmap_nproc"
    GMAP_NPROC_DEFAULT = 24


def add_map_isoforms_io_arguments(arg_parser):
    """Add io arguments"""
    helpstr = "Input isoforms in FASTA/FASTQ/Dataset format."
    arg_parser.add_argument("input_filename", type=str,
                            help=helpstr)

    helpstr = "Output SAM file mapping isoforms to reference using GMAP."
    arg_parser.add_argument("sam_filename", type=str,
                            help=helpstr)

    return arg_parser


def add_gmap_arguments(arg_parser):
    """Add gmap arguments"""
    helpstr = "GMAP DB name (default: %s)" % Constants.GMAP_NAME_DEFAULT
    arg_parser.add_argument("--gmap_name", type=str,
                            default=Constants.GMAP_NAME_DEFAULT,
                            help=helpstr)

    helpstr = "GMAP DB location (default: %s)" % Constants.GMAP_DB_DEFAULT
    arg_parser.add_argument("--gmap_db", type=str,
                            default=Constants.GMAP_DB_DEFAULT, help=helpstr)

    helpstr = "GMAP nproc (default: %s)" % Constants.GMAP_NPROC_DEFAULT
    arg_parser.add_argument("--gmap_nproc", type=int,
                            default=Constants.GMAP_NPROC_DEFAULT, help=helpstr)

    return arg_parser


def args_runner(args):
    """Run given input args"""
    try:
        map_isoforms_and_sort(input_filename=args.input_filename,
                              sam_filename=args.sam_filename,
                              gmap_db_dir=args.gmap_db,
                              gmap_db_name=args.gmap_name,
                              gmap_nproc=args.gmap_nproc)
    except Exception as e:
        raise RuntimeError("%s failed: %s" % (__file__, str(e)))
    return 0


def resolved_tool_contract_runner(rtc):
    """Run given a resolved tool contract"""
    try:
        map_isoforms_and_sort(input_filename=rtc.task.input_files[0],
                              sam_filename=rtc.task.output_files[0],
                              gmap_db_dir=rtc.task.options[Constants.GMAP_DB_ID],
                              gmap_db_name=rtc.task.options[Constants.GMAP_NAME_ID],
                              gmap_nproc=rtc.task.options[Constants.GMAP_NPROC_ID])
    except Exception as e:
        raise RuntimeError("%s failed: %s" % (__file__, str(e)))
    return 0


def get_contract_parser():
    """Get tool contract parser.
    Input:
        idx 0 - gmap input contigset
    Output:
        idx 0 - gmap output SAM
    """
    p = get_base_contract_parser(Constants, default_level="DEBUG")

    # argument parser
    add_map_isoforms_io_arguments(p.arg_parser.parser)
    add_gmap_arguments(p.arg_parser.parser)

    # tool contract parser
    tcp = p.tool_contract_parser
    tcp.add_input_file_type(FileTypes.DS_CONTIG, "gmap_input_ds", "ContigSet In",
                            "Gmap input ContigSet file") # input 0

    tcp.add_output_file_type(FileTypes.SAM, "gmap_output_sam",
                             name="SAM file", description="Gmap output sam",
                             default_name="gmap_output")

    tcp.add_str(option_id=Constants.GMAP_DB_ID, option_str="gmap_db",
                default=Constants.GMAP_DB_DEFAULT,
                name="GMAP DB Path", description="GMAP DB root directory")

    tcp.add_str(option_id=Constants.GMAP_NAME_ID, option_str="gmap_name",
                default=Constants.GMAP_NAME_DEFAULT,
                name="GMAP DB Name", description="GMAP DB Name")

    tcp.add_int(option_id=Constants.GMAP_NPROC_ID, option_str="gmap_nproc",
                default=Constants.GMAP_NPROC_DEFAULT,
                name="GMAP nproc", description="GMAP nproc")
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
