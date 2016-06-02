#!/usr/bin/env python

"""
Utils for mapping isoforms to reference genomes and sort.
"""

import os.path as op
import logging
from pbtranscript.Utils import execute, rmpath
from pbtranscript.io import ContigSetReaderWrapper

logger = logging.getLogger(op.basename(__file__))

def copy_sam_header(in_sam, out_sam):
    """Copy headers of input sam to output sam."""
    with open(in_sam, 'r') as reader, \
        open(out_sam, 'w') as writer:
        for line in reader:
            if line.startswith('@'):
                writer.write(line)
            else:
                break


def map_isoforms_and_sort(input_filename, sam_filename,
                          gmap_db_dir, gmap_db_name, gmap_nproc):
    """
    Map isoforms to references by gmap, generate a sam output and sort sam.
    Parameters:
        input_filename -- input isoforms. e.g., hq_isoforms.fasta|fastq|xml
        sam_filename -- output sam file, produced by gmap and sorted.
        gmap_db_dir -- gmap database directory
        gmap_db_name -- gmap database name
        gmap_nproc -- gmap nproc
    """
    unsorted_sam_filename = sam_filename + ".tmp"
    log_filename = sam_filename + ".log"

    gmap_input_filename = input_filename
    if input_filename.endswith('.xml'):
        # must consolidate dataset xml to FASTA/FASTQ
        w = ContigSetReaderWrapper(input_filename)
        gmap_input_filename = w.consolidate(out_prefix=sam_filename+'.input')
    if not op.exists(gmap_input_filename):
        raise IOError("Gmap input file %s does not exists" % gmap_input_filename)

    cmd_args = ['gmap', '-D {d}'.format(d=gmap_db_dir),
                '-d {name}'.format(name=gmap_db_name),
                '-t {nproc}'.format(nproc=gmap_nproc),
                '-n 0',
                '-z sense_force',
                '--cross-species',
                '-f samse',
                gmap_input_filename,
                '>', unsorted_sam_filename,
                '2>{log}'.format(log=log_filename)]
    # Call gmap to map isoforms to reference and output sam.
    execute(' '.join(cmd_args))

    # Copy SAM headers
    copy_sam_header(in_sam=unsorted_sam_filename,
                    out_sam=sam_filename)

    # Call sort to sort gmap output sam file
    cmd_args = ['sort', '-k 3,3', '-k 4,4n', unsorted_sam_filename,
                '| grep -v \'^@\'', '>>', sam_filename]

    execute(' '.join(cmd_args))

    # remove intermediate unsorted sam file.
    rmpath(unsorted_sam_filename)

