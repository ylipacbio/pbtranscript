#!/usr/bin/env python
"""
Specialized scatter task which spawns partial_chunks.pickle.
"""

import logging
import os.path as op
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes, PipelineChunk
from pbcommand.pb_io import write_pipeline_chunks

from pbtranscript.tasks.TPickles import ChunkTasksPickle, PartialChunkTask


log = logging.getLogger(__name__)


class Constants(object):
    """Constans used in pbtranscript.tasks.scatter_ice_partial_cluster_bins"""
    TOOL_ID = "pbtranscript.tasks.scatter_ice_partial_cluster_bins"
    DEFAULT_NCHUNKS = 24
    VERSION = "0.1.0"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID
    CHUNK_KEYS = ('$chunk.partial_chunk_pickle_id',
                  '$chunk.sentinel_id', '$chunk.ccs_id')


def get_contract_parser():
    """
    input:
      idx 0: partial_chunk_pickle_id
      idx 1: sentinel txt file
      idx 2: ccs_id
    output:
      idx 0: chunk json
    """
    p = get_scatter_pbparser(Constants.TOOL_ID, Constants.VERSION,
                             "Scatter Ice Partial Chunks",
                             __doc__, Constants.DRIVER_EXE,
                             chunk_keys=Constants.CHUNK_KEYS,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.PICKLE, "partial_chunk_pickle",
                          "PICKLE", "Partial Chunk Tasks Pickle") # input idx 0
    p.add_input_file_type(FileTypes.TXT, "partial_sentinel_in", "Sentinel In",
                          "Setinel file") # input idx 1
    p.add_input_file_type(FileTypes.DS_CCS, "ccs_in", "ConsensusReadSet In",
                          "PacBio ConsensusReadSet") # input idx 2
    p.add_output_file_type(FileTypes.CHUNK, "cjson_out",
                           "Chunk JSON Partial Tasks",
                           "Chunked JSON Partial Tasks",
                           "ice_partial.chunked")
    # max nchunks for this specific task
    p.add_int("pbsmrtpipe.task_options.dev_scatter_max_nchunks", "max_nchunks",
              Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")
    return p


def run_main(partial_chunks_pickle_file, sentinel_file,
             ccs_file, output_json_file, max_nchunks):
    """
    Spawn partial Chunk Tasks in pickle.
    Parameters:
      partial_chunks_pickle_file -- ChunkTasksPickle of PartialChunkTask objects
      ccs_file -- ccs dataset
      sentinel_file -- sentinel file to connect pbsmrtpipe tasks
      output_json -- chunk.json
    """
    p = ChunkTasksPickle.read(partial_chunks_pickle_file)
    assert all([isinstance(r, PartialChunkTask) for r in p])
    out_dir = op.dirname(output_json_file)

    # sort and group tasks
    groups = p.sort_and_group_tasks(max_nchunks=max_nchunks)

    # Writing chunk.json
    base_name = "spawned_partial_chunk"
    chunks = []
    spawned_pickles = []
    for group_index in range(0, len(groups)):
        chunk_id = "_".join([base_name, 'group', str(group_index)])
        spawned_pickle_file = op.join(out_dir, chunk_id + ".pickle")
        d = {Constants.CHUNK_KEYS[0]: spawned_pickle_file,
             Constants.CHUNK_KEYS[1]: sentinel_file,
             Constants.CHUNK_KEYS[2]: ccs_file}
        c = PipelineChunk(chunk_id, **d)
        chunks.append(c)
        spawned_pickles.append(spawned_pickle_file)

    log.info("Spawning %s into %d files", partial_chunks_pickle_file, len(groups))
    p.spawn_pickles_by_groups(groups=groups, out_pickle_fns=spawned_pickles)
    log.debug("Spawned files: %s.", ", ".join(spawned_pickles))

    log.info("Writing chunk.json to %s", output_json_file)
    write_pipeline_chunks(chunks, output_json_file,
                          "created by %s" % Constants.TOOL_ID)
    return 0


def args_run(args):
    """Args runner."""
    raise NotImplementedError()


def rtc_runner(rtc):
    """Resolved tool contract runner."""
    return run_main(partial_chunks_pickle_file=rtc.task.input_files[0],
                    sentinel_file=rtc.task.input_files[1],
                    ccs_file=rtc.task.input_files[2],
                    output_json_file=rtc.task.output_files[0],
                    max_nchunks=rtc.task.max_nchunks)


def main():
    """Main"""
    mp = get_contract_parser()
    return pbparser_runner(sys.argv[1:],
                           mp,
                           args_run,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
