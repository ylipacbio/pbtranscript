#!/usr/bin/env python

"""
Utils for combining output files from cluster bins.
"""

import os.path as op
import logging
import cPickle

from pbcore.io import FastqReader, FastqWriter, FastaWriter

from pbtranscript.io.ContigSetReaderWrapper import ContigSetReaderWrapper
from pbtranscript.Utils import mkdir


class CombinedFiles(object):
    """Combine cluster bin output files and save to combined_dir."""
    def __init__(self, combined_dir):
        self.combined_dir = combined_dir
        mkdir(self.combined_dir)

    @property
    def all_hq_fa(self):
        """Return combined polished lq isoforms in fasta"""
        return op.join(self.combined_dir, 'all.polished_hq.fasta')

    @property
    def all_hq_fq(self):
        """Return combined polished lq isoforms in fastq"""
        return op.join(self.combined_dir, 'all.polished_hq.fastq')

    @property
    def all_lq_fa(self):
        """Return combined polished lq isoforms in fasta"""
        return op.join(self.combined_dir, 'all.polished_lq.fasta')

    @property
    def all_lq_fq(self):
        """Return combined polished lq isoforms in fastq"""
        return op.join(self.combined_dir, 'all.polished_lq.fastq')

    @property
    def hq_lq_prefix_dict_pickle(self):
        """hq or lq prefix --> file path"""
        return op.join(self.combined_dir, 'all.hq_lq_pre_dict.pickle')

    @property
    def all_consensus_isoforms_fa(self):
        """Return combined unpolished consensus isoforms."""
        return op.join(self.combined_dir, 'all.consensus_isoforms.fasta')

    @property
    def all_cluster_report_fn(self):
        """Return combined cluster report."""
        return op.join(self.combined_dir, 'all.cluster_report.csv')

    @property
    def all_cluster_summary_fn(self):
        """Return combined cluster summary."""
        return op.join(self.combined_dir, 'all.cluster_summary.json')


def combined_prefix(cluster_bin_index, isoform_type, sample_name):
    """Return i{cluster_bin_index}_{isoform_type}_{sample_name}|"""
    assert isoform_type in ["HQ", "LQ", "ICE"]
    return "i{i}_{t}_{s}|".format(i=cluster_bin_index, t=isoform_type, s=sample_name)


def combined_cid_ice_name(name, cluster_bin_index, sample_name):
    """e.g., c1 --> i0_ICE_sample1|c1
    """
    return "{p}{n}".format(p=combined_prefix(cluster_bin_index=cluster_bin_index,
                                             isoform_type="ICE",
                                             sample_name=sample_name), n=name)

def combined_cid_hq_name(name, cluster_bin_index, sample_name):
    """e.g., c1 --> i0_HQ_sample1|c1
    """
    return "{p}{n}".format(p=combined_prefix(cluster_bin_index=cluster_bin_index,
                                             isoform_type="HQ",
                                             sample_name=sample_name), n=name)

def combined_cid_lq_name(name, cluster_bin_index, sample_name):
    """e.g., c1 --> i0_LQ_sample1|c1
    """
    return "{p}{n}".format(p=combined_prefix(cluster_bin_index=cluster_bin_index,
                                             isoform_type="LQ",
                                             sample_name=sample_name), n=name)


def combine_polished_isoforms(split_indices, split_hq_fns, split_lq_fns,
                              combined_hq_fa, combined_hq_fq,
                              combined_lq_fa, combined_lq_fq,
                              hq_lq_prefix_dict_pickle, sample_name):
    """Combine split hq (lq) files and save to combined_dir.
    Dumping hq|lq prefix dictionary to pickle.
    Return an instance of CombinedFiles.
    Parameters:
      split_indices -- indices of splitted cluster bins.
      split_hq_fns -- hq files, #['*/all_quivered_hq.100_30_0.99.fastq', ...]
      split_lq_fns -- lq files, #['all_quivered_lq.fastq', ...]
    """
    assert len(split_indices) == len(split_hq_fns)
    assert len(split_indices) == len(split_lq_fns)
    assert all([f.endswith(".fastq") for f in split_hq_fns + split_lq_fns])

    hq_pre_dict, lq_pre_dict = {}, {}

    hq_fa_writer = FastaWriter(combined_hq_fa)
    hq_fq_writer = FastqWriter(combined_hq_fq)
    lq_fa_writer = FastaWriter(combined_lq_fa)
    lq_fq_writer = FastqWriter(combined_lq_fq)

    for i, split_hq, split_lq in zip(split_indices, split_hq_fns, split_lq_fns):
        logging.debug("Adding prefix i%s_| to %s, %s", str(i), split_hq, split_lq)
        hq_prefix = combined_prefix(cluster_bin_index=i, isoform_type="HQ",
                                    sample_name=sample_name)
        lq_prefix = combined_prefix(cluster_bin_index=i, isoform_type="LQ",
                                    sample_name=sample_name)

        hq_pre_dict[hq_prefix] = op.dirname(op.abspath(split_hq))
        lq_pre_dict[lq_prefix] = op.dirname(op.abspath(split_lq))

        with FastqReader(split_hq) as reader:
            for read in reader:
                name = combined_cid_hq_name(cluster_bin_index=i,
                                            name=read.name, sample_name=sample_name)
                hq_fa_writer.writeRecord(name, read.sequence[:])
                hq_fq_writer.writeRecord(name, read.sequence[:], read.quality)

        with FastqReader(split_lq) as reader:
            for read in reader:
                name = combined_cid_lq_name(cluster_bin_index=i,
                                            name=read.name, sample_name=sample_name)
                lq_fa_writer.writeRecord(name, read.sequence[:])
                lq_fq_writer.writeRecord(name, read.sequence[:], read.quality)
    hq_fa_writer.close()
    hq_fq_writer.close()
    lq_fa_writer.close()
    lq_fq_writer.close()
    logging.info("HQ polished output combined to:%s", combined_hq_fq)
    logging.info("LQ polished output combined to:%s", combined_lq_fq)

    logging.info("Dumping hq|lq prefix dictionary to:%s", hq_lq_prefix_dict_pickle)
    with open(hq_lq_prefix_dict_pickle, 'wb') as writer:
        cPickle.dump({'HQ': hq_pre_dict, 'LQ': lq_pre_dict}, writer)


def combine_consensus_isoforms(split_indices, split_files,
                               combined_consensus_isoforms_fa,
                               sample_name):
    """
    Parameters:
      split_indices -- indices of splitted cluster bins.
      split_files -- consensus isoforms in each splitted cluster bin.
    """
    assert len(split_indices) == len(split_files)
    writer = FastaWriter(combined_consensus_isoforms_fa)
    for i, split_fn in zip(split_indices, split_files):
        logging.debug("Adding prefix i%s to %s.", str(i), split_fn)
        with ContigSetReaderWrapper(split_fn) as reader:
            for read in reader:
                name = combined_cid_ice_name(name=read.name, cluster_bin_index=i,
                                             sample_name=sample_name)
                writer.writeRecord(name, read.sequence[:])
    writer.close()
    logging.info("Consensus isoforms output combined to:%s",
                 combined_consensus_isoforms_fa)


def write_combined_cluster_report(split_indices, split_uc_pickles,
                                  split_partial_uc_pickles, report_fn,
                                  sample_name):
    """
    Write a CSV report to report_fn, each line contains three columns:
        cluster_id, read_id and read_type
    e.g., i0_ICE_samplename|c1 m12345/123/0_1000 FL

    Parameters:
      split_indices -- indices of splitted cluster bins.
      split_uc_pickles -- uc pickle (output/final.pickle) in
                          each splitted cluster bin.
      split_partial_uc_pickles -- partial uc pickle
                          (output/map_noFL/nfl.all.partial.pickle)
                          in each splitted cluster bin.
    """
    assert len(split_indices) == len(split_uc_pickles)
    assert len(split_indices) == len(split_partial_uc_pickles)

    with open(report_fn, 'w') as f:
        f.write("cluster_id,read_id,read_type\n")
        for i, uc_pickle, partial_uc_pickle in zip(split_indices,
                                                   split_uc_pickles,
                                                   split_partial_uc_pickles):
            logging.info("Processing uc pickle %s and partial uc pickle %s",
                         uc_pickle, partial_uc_pickle)
            uc = cPickle.load(open(uc_pickle, 'rb'))['uc']
            partial_uc = cPickle.load(open(partial_uc_pickle, 'rb'))['partial_uc']
            for c in uc.keys():
                for r in uc[c]:
                    cid = combined_cid_ice_name(name="c{c}".format(c=c),
                                                cluster_bin_index=i,
                                                sample_name=sample_name)
                    f.write("{cid},{r},FL\n".format(cid=cid, r=r))
                if partial_uc is not None and c in partial_uc.keys():
                    for r in partial_uc[c]:
                        f.write("{cid},{r},NonFL\n".format(cid=cid, r=r))
