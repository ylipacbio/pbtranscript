#!/usr/bin/env python
"""
Filter collapsed isoforms by count or subset.
"""

from pbcore.io import FastaReader, FastaWriter, FastqReader, FastqWriter
from pbtranscript.io import GroupReader, AbundanceReader, AbundanceWriter, \
        CollapseGffReader, CollapseGffWriter, SampleIsoformName, parse_ds_filename


__author__ = 'etseng@pacb.com'


__all__ = ["filter_by_count"]


def filter_by_count(in_group_filename,
                    in_abundance_filename, in_gff_filename, in_rep_filename,
                    out_abundance_filename, out_gff_filename, out_rep_filename,
                    min_count):
    """Filter input (collapsed) isoforms in in_rep_filename by count.
    Parameters:
      in_group_filename -- collapsed isoforms' pbid --> clusters
      in_abundance_filename -- collapsed isoforms' pbid, count_fl, count_nfl, ...
      in_gff_filename -- collapsed isoforms' pbid, chr, strand, start, end, ...
      in_rep_filename -- representative sequences of collapsed isoforms
      min_count -- min value of count_fl and count_nfl to classify a collapsed isoform as good
    """
    in_suffix = parse_ds_filename(in_rep_filename)[1]
    out_suffix = parse_ds_filename(out_rep_filename)[1]
    if in_suffix != out_suffix:
        raise ValueError("Format of input %s and output %s must match." %
                         (in_rep_filename, out_rep_filename))
    if in_suffix not in ("fasta", "fastq"):
        raise ValueError("Format of input %s and output %s must be either FASTA or FASTQ." %
                         (in_rep_filename, out_rep_filename))

    # read group
    group_max_count_fl = {}
    group_max_count_nfl = {}
    with GroupReader(in_group_filename) as g_reader:
        for g in g_reader:
            pbid, members = g.name, g.members
            group_max_count_fl[pbid] = 0
            group_max_count_nfl[pbid] = 0
            for m in members:
                s = SampleIsoformName.fromString(m)
                group_max_count_fl[pbid] = max(group_max_count_fl[pbid], s.num_fl)
                group_max_count_nfl[pbid] = max(group_max_count_nfl[pbid], s.num_nfl)

    # read abundance to decide good collapsed isoforms based on count
    good = [r.pbid for r in AbundanceReader(in_abundance_filename)
            if r.count_fl >= min_count and group_max_count_fl[r.pbid] >= min_count]

    # then read gff, and write good gff record.
    with CollapseGffWriter(out_gff_filename) as gff_writer:
        for r in CollapseGffReader(in_gff_filename):
            if r.seqid in good:
                gff_writer.writeRecord(r)

    # next read fasta/fastq, and write good fasta/fastq record.
    rep_reader = FastaReader(in_rep_filename) if in_suffix == "fasta" \
                 else FastqReader(in_rep_filename)
    rep_writer = FastaWriter(out_rep_filename) if in_suffix == "fasta" \
                 else FastqWriter(out_rep_filename)
    for r in rep_reader:
        # r.name e.g., PB.1.1|PB.1.1:10712-11643(+)|i0_HQ_sample18ba5d|c1543/f8p1/465
        if r.name.split('|')[0] in good:
            rep_writer.writeRecord(r)

    # finally write filtered abundance report.
    with AbundanceReader(in_abundance_filename) as a_reader, \
        AbundanceWriter(out_abundance_filename, comments=a_reader.comments) as a_writer:
        for r in a_reader:
            if r.pbid in good:
                a_writer.writeRecord(r)
