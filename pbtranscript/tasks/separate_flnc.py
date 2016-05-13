# wraps 'pbtranscript.tasks.separate_flnc'

"""
Separate full-length-non-chimeric ccs reads into multiple
bins based on either primer or read length.
Inputs: isoseq_flnc.fasta|contigset.xml
Outputs:
        <out_pickle>
        <root_dir>/<bin>/isoseq_flnc.fasta|contigset.xml
"""

import logging
import os.path as op
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptOptions import get_base_contract_parser
from pbtranscript.separate_flnc import Constants, \
        add_separate_flnc_arguments, args_runner, get_arg_parser

log = logging.getLogger(__name__)


def get_contract_parser():
    """Return a tool contract parser for separate_flnc."""
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    add_separate_flnc_arguments(p)
    return p


def resolved_tool_contract_to_args(resolved_tool_contract):
    """
    Convert an input resolved tool contract to args.
    """
    rtc = resolved_tool_contract
    # rtc has only one input, idx 0, isoseq_flnc.contigset.xml
    # rtc has only one output, idx 0, separate_flnc.pickle
    args = [rtc.task.input_files[0],
            op.dirname(rtc.task.output_files[0]),
            "--out_pickle", rtc.task.output_files[0]]
    by_primer = rtc.task.options[Constants.BIN_BY_PRIMER_ID] is True
    if by_primer: # by primer, over write by_manual, by_size_kb
        args.append("--bin_by_primer")
    else:
        bin_manual = rtc.task.options[Constants.BIN_MANUAL_ID]
        if (isinstance(bin_manual, str) or isinstance(bin_manual, unicode)) \
            and len(str(bin_manual).translate(None, '[]() ')) > 0:
            log.info("Setting bin manually, bin_manual is %s" % bin_manual)
            args.extend(["--bin_manual", bin_manual])
        else: # bin_manual is None
            args.extend(["--bin_size_kb", str(rtc.task.options[Constants.BIN_SIZE_KB_ID])])
    return get_arg_parser().parse_args(args)


def resolved_tool_contract_runner(resolved_tool_contract):
    """Given an input resolved tool contract, convert it to args and
       call args_runner."""
    args = resolved_tool_contract_to_args(resolved_tool_contract)
    return args_runner(args)


def main(argv=sys.argv[1:]):
    """Main"""
    return pbparser_runner(
        argv=argv,
        parser=get_contract_parser(),
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
