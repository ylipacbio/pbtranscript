#!/usr/bin/env python

"""Access data files in SMRTLink IsoSeq job"""

import logging
import os.path as op
import json
from pbcore.io import ContigSet, FastaWriter, FastqWriter

log = logging.getLogger(__file__)

__all__ = [
    "smrtlink_dir",
    "consolidate_xml",
    "json_to_attr_dict"
]

def smrtlink_dir(smrtlink_host, job_id):
    """Given smrtlink host and job id, return path to a smrtlink job dir."""
    host_jobs_root = ""
    if smrtlink_host == "alpha":
        host_jobs_root = "/pbi/dept/secondary/siv/smrtlink/smrtlink-alpha/jobs-root/"
    elif smrtlink_host == "beta":
        host_jobs_root = "/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/jobs-root/"
    else:
        raise ValueError("SMRTLink host %s is not supported." % smrtlink_host)

    return op.join(host_jobs_root,
                   "{0:03d}".format(int(int(job_id)/1000)),
                   "{0:06d}".format(int(job_id)))


def consolidate_xml(src, dst):
    """Convert input dataset to output fasta|fastq"""
    w = None
    if dst.endswith(".fa") or dst.endswith(".fasta"):
        w = FastaWriter(dst)
        for r in ContigSet(src):
            w.writeRecord(r)
        w.close()
    elif dst.endswith(".fq") or dst.endswith(".fastq"):
        w = FastqWriter(dst)
        for r in ContigSet(src):
            w.writeRecord(r)
        w.close()
    else:
        raise ValueError("Output file %s must be either fasta or fastq", dst)


def load_json(json_fn):
    """"Read a json to a dict"""
    with open(json_fn) as json_data:
        return json.load(json_data)


def json_to_attr_dict(json_fn):
    """Convert a json report to a dict {id: value}."""
    ret = dict()
    d = load_json(json_fn)
    try:
        for attr in d['attributes']:
            ret[attr['id']] = attr['value']
    except ValueError:
        log.warning("Warning: could not get attributes from json report %s", json_fn)
    return ret
