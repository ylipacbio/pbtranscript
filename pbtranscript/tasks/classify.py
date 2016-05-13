
# wraps 'pbtranscript classify'

"""
Classifies reads from a fasta/q file.  For each read, identify whether it is
full length, whether 5', 3' and poly A tail have been found. The input is a
ConsensusRead dataset.
"""

import logging
import os.path as op
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptRunner import PBTranscript
from pbtranscript.PBTranscriptOptions import (BaseConstants,
                                              get_base_contract_parser,
                                              get_argument_parser,
                                              add_classify_arguments)

log = logging.getLogger(__name__)

class Constants(BaseConstants):
    TOOL_ID = "pbtranscript.tasks.classify"
    DRIVER_EXE = "python -m pbtranscript.tasks.classify --resolved-tool-contract"
    PARSER_DESC = __doc__


def get_contract_parser():
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    add_classify_arguments(p)
    return p


def args_runner(args):
    return PBTranscript(args, subCommand="classify").start()


def resolved_tool_contract_to_args(resolved_tool_contract):
    rtc = resolved_tool_contract
    args = [
        "--verbose",
        "classify",
        resolved_tool_contract.task.input_files[0],
        resolved_tool_contract.task.output_files[0],
        "--flnc", resolved_tool_contract.task.output_files[1],
        "--nfl", resolved_tool_contract.task.output_files[2],
        "--summary", resolved_tool_contract.task.output_files[3],  # JSON
        "--report", resolved_tool_contract.task.output_files[4],  # CSV
        "--min_seq_len", str(rtc.task.options[Constants.MIN_SEQ_LEN_ID]),
        "--cpus", str(resolved_tool_contract.task.nproc),
        "--outDir", op.dirname(rtc.task.output_files[0]),
        "--ignore-empty-output",
    ]
    if rtc.task.options[Constants.IGNORE_POLYA_ID]:
        args.append("--ignore_polyA")
    return get_argument_parser().parse_args(args)


def resolved_tool_contract_runner(resolved_tool_contract):
    args = resolved_tool_contract_to_args(resolved_tool_contract)
    return args_runner(args)


def main(argv=sys.argv[1:]):
    return pbparser_runner(
        argv=argv,
        parser=get_contract_parser(),
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
