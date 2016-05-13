
import tempfile
import unittest
import logging
import cPickle
import re
import os.path as op
import os

from pbcoretools.chunking.chunk_utils import write_chunks_to_json
from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.pb_io.report import load_report_from_json
from pbcommand.models import PipelineChunk
import pbcommand.testkit.core
from pbcore.io import ContigSet

TEST_DIR = op.dirname(op.dirname(__file__))
DATA = op.join(TEST_DIR, "data", "report_gather")


class TestReportGather(pbcommand.testkit.core.PbTestGatherApp):
    DRIVER_BASE = "python -m pbtranscript.tasks.gather_isoseq_classify_report"
    CHUNK_KEY = "$chunk.report_id"
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".chunks.json").name
    ]

    def setUp(self):
        data_files = [op.join(DATA, fn) for fn in os.listdir(DATA)
                      if fn.startswith("summary")]
        chunks = [PipelineChunk(chunk_id="chunk_data_{i}".format(i=i),
                                **({self.CHUNK_KEY:fn}))
                  for i, fn in enumerate(data_files)]
        write_chunks_to_json(chunks, self.INPUT_FILES[0])

    def run_after(self, rtc, output_dir):
        rpt = load_report_from_json(rtc.task.output_files[0])
        for a in rpt.attributes:
            if a.id == "avg_flnc_len":
                self.assertEqual(a.value, 1142)
                break
        else:
            self.fail("avg_flnc_len not found")
