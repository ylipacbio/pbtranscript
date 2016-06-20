#!/usr/bin/env python

"""
Map isoforms to reference genomes and sort.
"""
import sys
import logging

from pbcore.io import GmapReferenceSet

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptOptions import get_base_contract_parser
from pbtranscript.collapsing.CollapsingUtils import map_isoforms_and_sort

log = logging.getLogger(__name__)

class Constants(object):
    """Constants used in tool contract."""
    TOOL_ID = "pbtranscript.tasks.map_isoforms_to_genome"
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

    helpstr = "GMAP Reference Set file to overwrite GMAP DB name and location."
    arg_parser.add_argument("--gmap_ds", type=str, default=None, help=helpstr)

    helpstr = "GMAP nproc (default: %s)" % Constants.GMAP_NPROC_DEFAULT
    arg_parser.add_argument("--gmap_nproc", type=int,
                            default=Constants.GMAP_NPROC_DEFAULT, help=helpstr)

    return arg_parser


def gmap_db_and_name_from_ds(gmap_ds_filename):
    """Return gmap db dir and gmap db name"""
    gmap_ds = GmapReferenceSet(gmap_ds_filename)
    return (gmap_ds.gmap.resourceId, gmap_ds.gmap.name)


def args_runner(args):
    """Run given input args.
    e.g.,
    map_isoforms_to_genome.py hq_isoforms.fastq out.sam --gmap_db=<path-to-db> --gmap_name=<name>
    map_isoforms_to_genome.py hq_isoforms.fastq out.sam --gmap_ds=<path-to-xml>
    """
    gmap_db_dir, gmap_db_name = args.gmap_db, args.gmap_name
    if args.gmap_ds is not None:
        gmap_db_dir, gmap_db_name = gmap_db_and_name_from_ds(args.gmap_ds)

    map_isoforms_and_sort(input_filename=args.input_filename,
                          sam_filename=args.sam_filename,
                          gmap_db_dir=gmap_db_dir,
                          gmap_db_name=gmap_db_name,
                          gmap_nproc=args.gmap_nproc)
    return 0


def resolved_tool_contract_runner(rtc):
    """Run given a resolved tool contract"""
    gmap_db_dir, gmap_db_name = gmap_db_and_name_from_ds(rtc.task.input_files[1])
    map_isoforms_and_sort(input_filename=rtc.task.input_files[0],
                          sam_filename=rtc.task.output_files[0],
                          gmap_db_dir=gmap_db_dir,
                          gmap_db_name=gmap_db_name,
                          gmap_nproc=rtc.task.options[Constants.GMAP_NPROC_ID])
    return 0


def get_contract_parser():
    """Get tool contract parser.
    Input:
        idx 0 - HQ isoforms fastq as gmap input
        idx 1 - gmap reference set
    Output:
        idx 0 - gmap output SAM
    """
    p = get_base_contract_parser(Constants, default_level="DEBUG")

    # argument parser
    add_map_isoforms_io_arguments(p.arg_parser.parser)
    add_gmap_arguments(p.arg_parser.parser)

    # tool contract parser
    tcp = p.tool_contract_parser
    tcp.add_input_file_type(FileTypes.FASTQ, "hq_isoforms_fastq", "FASTQ In",
                            "HQ isoforms FASTQ file") # input 0
    tcp.add_input_file_type(FileTypes.DS_GMAP_REF, "gmap_referenceset", "GmapReferenceSet In",
                            "Gmap reference set file") # input 1
    tcp.add_output_file_type(FileTypes.SAM, "gmap_output_sam",
                             name="SAM file", description="Gmap output sam",
                             default_name="gmap_output")

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
