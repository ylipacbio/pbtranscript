#!/usr/bin/env python
"""
MegaPBTree, non-redundant set of gene annotations for
combining collapsed GFFs from mulitiple smaples.

The original file name: combine_abundance_across_samples.py
"""
import copy
from collections import defaultdict
from pbtranscript.collapsing import IntervalTree, compare_fuzzy_junctions
from pbtranscript.collapsing.cluster import ClusterTree
from pbtranscript.io import (CollapseGffReader, CollapseGffWriter, CollapseGffRecord,
                             GmapRecord, GroupReader, GroupWriter, GroupRecord,
                             MegaInfoWriter, MergeGroupOperation)

__author__ = 'etseng@pacificbiosciences.com'

__all__ = ['MegaPBTree']


def read_gff_as_interval_tree(gff_filename):
    """
    Read a collapsed GFF file into an IntervalTree
    """
    tree = defaultdict(
        lambda: {'+': IntervalTree(), '-': IntervalTree()})  # chr --> strand --> tree
    for r in CollapseGffReader(gff_filename):
        tree[r.chr][r.strand].insert(r.start, r.end, r)
    return tree


class MegaPBTree(object):
    """
    Structure for maintaining a non-redundant set of gene annotations
    Used to combine with different collapsed GFFs from different samples
    """

    def __init__(self, gff_filename, group_filename, self_prefix=None, max_fuzzy_junction=0):
        self.gff_filename = gff_filename
        self.group_filename = group_filename
        self.self_prefix = self_prefix
        self.max_fuzzy_junction = max_fuzzy_junction

        self.record_d = dict((r.seqid, r)
                             for r in CollapseGffReader(gff_filename))
        self.tree = read_gff_as_interval_tree(
            gff_filename=self.gff_filename)  # chr --> strand -->tree
        # ex: PB.1.1 --> [ RatHeart|i3_c123.... ]
        self.group_info = MegaPBTree.read_group(
            self.group_filename, self.self_prefix)

        # keep track of gff|group files that has been added.
        self._sample_prefixes = []
        self._group_filenames = []
        self._gff_filenames = []
        self._add_sample_files(
            gff_filename=gff_filename, group_filename=group_filename, sample_prefix="first_sample")

    def __str__(self):
        ret = ["MegaPBTree of %s samples" % len(self.gff_filename)]
        for sample, gff_fn, group_fn in zip(self._sample_prefixes, self._gff_filenames, self._group_filenames):
            ret.extend(["sample %s\n" % sample, "\t%s",
                        gff_fn, "\t%s" % group_fn])
        return '\n'.join(ret)

    @staticmethod
    def read_group(group_filename, group_prefix):
        """read a group file and group_prefix to a dict
        if group_prefix is None: return {group.pbid --> group.members}
        else: {group.pbid --> [group_prefix+'|'+m for m in group.members]
        """
        return {group.name: group.members
                for group in GroupReader(group_filename, group_prefix)}

    def match_record_to_tree(self, r):
        """
        r --- GmapRecord
        tree --- dict of chromosome --> strand --> IntervalTree

        If exact match (every exon junction), return the matching GmapRecord
        Otherwise return None
        *NOTE*: the tree should be non-redundant so can return as soon as exact match is found!
        """
        assert isinstance(r, GmapRecord)
        matches = self.tree[r.chr][r.strand].find(r.start, r.end)
        for r2 in matches:
            #r.segments = r.ref_exons
            #r2.segments = r2.ref_exons
            # is a match!
            if compare_fuzzy_junctions(r.ref_exons, r2.ref_exons, self.max_fuzzy_junction) == 'exact':
                return r2
        return None

    def add_sample(self, gff_filename, group_filename, sample_prefix, o_gff_fn, o_group_fn, o_mega_fn):
        """Add one more sample to this MagaPBTree object.
        Read gff file to get collapsed isoforms from new sample,
        combine with existing collapsed isoforms and update tree.
        """
        self._add_sample_files(
            gff_filename=gff_filename, group_filename=group_filename, sample_prefix=sample_prefix)

        # list of (r1 if r2 is None | r2 if r1 is None | longer of r1 or r2 if
        # both not None)
        combined = []
        unmatched_recs = self.record_d.keys()

        for r in CollapseGffReader(gff_filename):
            match_rec = self.match_record_to_tree(r)
            if match_rec is not None:  # found a match! put longer of r1/r2 in
                combined.append((match_rec, r))
                try:
                    unmatched_recs.remove(match_rec.seqid)
                except ValueError:
                    pass  # already deleted, OK, this happens for single-exon transcripts
            else:  # r is not present in current tree
                combined.append((None, r))
        # put whatever is left from the tree in
        for seqid in unmatched_recs:
            combined.append((self.record_d[seqid], None))

        # create a ClusterTree to re-calc the loci/transcripts
        final_tree = defaultdict(
            lambda: {'+': ClusterTree(0, 0), '-': ClusterTree(0, 0)})
        for i, (r1, r2) in enumerate(combined):
            if r2 is None or (r1 is not None and r1.end - r1.start > r2.end - r2.start):
                final_tree[r1.chr][r1.strand].insert(r1.start, r1.end, i)
            else:
                final_tree[r2.chr][r2.strand].insert(r2.start, r2.end, i)

        self.write_cluster_tree_as_gff(
            final_tree, combined, group_filename, sample_prefix, o_gff_fn, o_group_fn, o_mega_fn)

    def write_cluster_tree_as_gff(self, cluster_tree, rec_list, group_filename2, sample_prefix2, o_gff_fn, o_group_fn, o_mega_fn):
        """
        Write ClusterTree (chr --> dict --> (start, end, rec_list_index)) as collapsedGFF format
        Returns --- a new group_info!!!
        rec_list --- a list of (r1, r2) where r1 and r2 are GmapRecord
        """
        group_info2 = MegaPBTree.read_group(group_filename2, sample_prefix2)
        new_group_info = {}

        gff_writer = open(o_gff_fn, 'w')
        group_writer = GroupWriter(o_group_fn)
        f_mgroup_writer = MegaInfoWriter(o_mega_fn, self.self_prefix, sample_prefix2)

        loci_index = 0
        chroms = cluster_tree.keys()
        chroms.sort()
        for k in chroms:
            for strand in ('+', '-'):
                for dummy_s, dummy_e, rec_indices in cluster_tree[k][strand].getregions():
                    loci_index += 1
                    isoform_index = 0
                    for i in rec_indices:
                        isoform_index += 1
                        gene_id = "PB.{i}".format(i=loci_index)
                        tID = "{gene_id}.{j}".format(
                            gene_id=gene_id, j=isoform_index)
                        r1, r2 = rec_list[i]
                        assert isinstance(r1, GmapRecord) or r1 is None
                        assert isinstance(r2, GmapRecord) or r2 is None
                        if r1 is None:  # r2 is not None
                            r = r2
                            new_group_info[tID] = group_info2[r2.seqid]
                        elif r2 is None:  # r1 is not None
                            r = r1
                            new_group_info[tID] = self.group_info[r1.seqid]
                        else:  # both r1, r2 are not empty
                            r = r1 if (r1.end - r1.start >
                                       r2.end - r2.start) else r2
                            new_group_info[tID] = self.group_info[r1.seqid] + \
                                group_info2[r2.seqid]

                        # write merged new group
                        group_writer.writeRecord(GroupRecord(name=tID, members=new_group_info[tID]))
                        # write group merge operation
                        f_mgroup_writer.writeRecord(MergeGroupOperation(pbid=tID, group1=r1, group2=r2))

                        gff_writer.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{tID}\";\n".format(
                            chr=k, s=r.start + 1, e=r.end, strand=strand, gene_id=gene_id, tID=tID))
                        for exon in r.ref_exons:
                            gff_writer.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{tID}\";\n".format(
                                chr=k, s=exon.start + 1, e=exon.end, strand=strand, gene_id=gene_id, tID=tID))

        gff_writer.close()
        group_writer.close()
        f_mgroup_writer.close()
        return new_group_info

    def _add_sample_files(self, gff_filename, group_filename, sample_prefix):
        """Keep track of gff|group files that has been added."""
        # keep track of gff|group filenames
        self._gff_filenames.append(gff_filename)
        self._group_filenames.append(group_filename)
        self._sample_prefixes.append(sample_prefix)
