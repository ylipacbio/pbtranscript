"""Test pbtranscript.collapsing.Branch."""
import unittest
import os.path as op
import cPickle
import filecmp
import numpy as np
from pbcore.io.GffIO import Gff3Record
from pbtranscript.Utils import rmpath, mkdir
from pbtranscript.io import ContigSetReaderWrapper, iter_gmap_sam, GroupWriter, CollapseGffWriter
from pbtranscript.collapsing import Branch, ContiVec, transfrag_to_contig, \
        exons_match_sam_record, compare_exon_matrix, get_fl_from_id, collapse_sam_records
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR, SIV_STD_DIR

_OUT_DIR_ = op.join(OUT_DIR, "test_branch")

SORTED_GMAP_SAM = op.join(SIV_DATA_DIR, 'test_branch', 'sorted-gmap-output.sam')
READS_DS = op.join(SIV_DATA_DIR, 'test_collapsing', 'gmap-input.fastq.contigset.xml')

def _get_sam_groups(ignored_ids_writer=None):
    """Returns grouped sam records read from SORTED_GMAP_SAM and READS_DS."""
    query_len_dict = ContigSetReaderWrapper.name_to_len_dict(READS_DS)
    groups = [g for g in iter_gmap_sam(sam_filename=SORTED_GMAP_SAM,
                                       query_len_dict=query_len_dict,
                                       min_aln_coverage=0.99, min_aln_identity=0.85,
                                       ignored_ids_writer=ignored_ids_writer)]
    return groups


def _get_contiVec_and_offset():
    """Returns contiVec and offset of groups[0]["+"]."""
    contivec_pickle_fn = op.join(SIV_DATA_DIR, 'test_branch', 'contiVec.pickle')
    a = cPickle.load(open(contivec_pickle_fn, 'rb'))
    return a['contiVec'], a['offset']


def _get_exons():
    """Returns exons of groups[0]["+"]."""
    contiVec, offset = _get_contiVec_and_offset()
    exons = contiVec.to_exons(offset=offset)
    return exons


class TEST_Branch(unittest.TestCase):
    """Test functions of pbtranscript.collapsing.Branch."""
    def setUp(self):
        """Define input and output file."""
        rmpath(_OUT_DIR_)
        mkdir(_OUT_DIR_)

    def test_transfrag_to_contig(self):
        """Test transfrag_to_contig, which takes a group of overlapping sam
        records as input and return (contiVec, offset, chrom, strand),
        where contiVec is a nparray vector of contig"""
        expected_contiVec, expected_offset = _get_contiVec_and_offset()
        groups = _get_sam_groups()

        contiVec, offset, chrom, strand = \
        transfrag_to_contig(groups[0]["+"], skip_5_exon_alt=True)
        self.assertEqual(offset, 10710)
        self.assertEqual(chrom, "SIRV1")
        self.assertEqual(strand, "+")
        self.assertEqual(contiVec, expected_contiVec)
        self.assertEqual(expected_offset, 10710)

    def test_ContiVec_to_exons(self):
        """Test ContiVec.to_exons"""
        exons = _get_exons() # contains 10 intervals
        p = []
        exons.traverse(p.append)
        expected_tree_0 = [(10710, 10712, 0), (10712, 10713, 1),
                           (10713, 10715, 2), (10715, 10791, 3),
                           (10791, 10791, 4), (10882, 11057, 5),
                           (11057, 11057, 6), (11434, 11638, 7),
                           (11638, 11640, 8), (11640, 11641, 9)]
        self.assertEqual([(node.start, node.end, node.interval.value) for node in p],
                         expected_tree_0)

    def test_exons_match_sam_record(self):
        """Test exons_match_sam_record, which takes a GMAP sam reocord and an exon tree
        (type IntervalUniqueTree, created by contiVec.to_exons) as input and
        returns a list of nodes this GMAP sam record corresponds to."""
        exons = _get_exons() # contains 10 intervals
        records = _get_sam_groups()[0]["+"] # contains 10 sam records
        self.assertEqual(len(records), 10)

        stuffs = [exons_match_sam_record(record=record, exons=exons) for record in records]
        # The very first sam record contains exons 0, 1, 2, 3, 5, 7, 8, 9, but not 4, 6.
        self.assertEqual([node.value for node in stuffs[0]], [0, 1, 2, 3, 5, 7, 8, 9])

        # The second sam record contains exons 1, 2, 3, 5, 7, 8, 9, but not 1, 4, 6.
        self.assertEqual([node.value for node in stuffs[1]], [1, 2, 3, 5, 7, 8, 9])

        # The third sam record contains exons 1, 2, 3, 5, 7, not 0, 4, 6, 8, 9
        self.assertEqual([node.value for node in stuffs[2]], [1, 2, 3, 5, 7])
        self.assertEqual([node.value for node in stuffs[3]], [1, 2, 3, 5, 7, 8])
        self.assertEqual([node.value for node in stuffs[4]], [1, 2, 3, 5, 7, 8, 9])
        self.assertEqual([node.value for node in stuffs[5]], [1, 2, 3, 5, 7, 8, 9])
        self.assertEqual([node.value for node in stuffs[6]], [1, 2, 3, 5, 7])
        self.assertEqual([node.value for node in stuffs[7]], [2, 3, 5, 7, 8, 9])
        self.assertEqual([node.value for node in stuffs[8]], [2, 3, 5, 7, 8])
        self.assertEqual([node.value for node in stuffs[9]], [3, 5, 7, 8, 9])

    def test_compare_exon_matrix(self):
        """
        test compare_exon_matrix, which takes two exon matrix (m1 and m2) as input and
        returns True if m1 and m2 can be merged and False otherwise.
        An exon matrix m is a 2-d array where m[0, i] is 1 if the 1-th exon will be used.
        """
        exons = _get_exons() # contains 10 intervals
        records = _get_sam_groups()[0]["+"] # contains 10 sam records
        stuffs = [exons_match_sam_record(record=record, exons=exons) for record in records]

        p = []
        exons.traverse(p.append)
        node_d = dict((x.interval.value, x) for x in p) # exon index --> exon node

        exon_all_indices = range(0, len(p))
        ms = [np.asarray([[1 if exon_index in [node.value for node in stuffs[record_index]] else 0
                           for exon_index in exon_all_indices]])
              for record_index in range(0, len(records))]

        self.assertTrue(np.all(ms[0] == np.asarray([[1, 1, 1, 1, 0, 1, 0, 1, 1, 1]])))
        self.assertTrue(np.all(ms[1] == np.asarray([[0, 1, 1, 1, 0, 1, 0, 1, 1, 1]])))
        self.assertTrue(np.all(ms[2] == np.asarray([[0, 1, 1, 1, 0, 1, 0, 1, 0, 0]])))
        self.assertTrue(np.all(ms[3] == np.asarray([[0, 1, 1, 1, 0, 1, 0, 1, 1, 0]])))
        self.assertTrue(np.all(ms[4] == np.asarray([[0, 1, 1, 1, 0, 1, 0, 1, 1, 1]])))
        self.assertTrue(np.all(ms[5] == np.asarray([[0, 1, 1, 1, 0, 1, 0, 1, 1, 1]])))
        self.assertTrue(np.all(ms[6] == np.asarray([[0, 1, 1, 1, 0, 1, 0, 1, 0, 0]])))
        self.assertTrue(np.all(ms[7] == np.asarray([[0, 0, 1, 1, 0, 1, 0, 1, 1, 1]])))
        self.assertTrue(np.all(ms[8] == np.asarray([[0, 0, 1, 1, 0, 1, 0, 1, 1, 0]])))
        self.assertTrue(np.all(ms[9] == np.asarray([[0, 0, 0, 1, 0, 1, 0, 1, 1, 1]])))

        for i in xrange(0, 10):
            for j in xrange(i, 10):
                self.assertTrue(compare_exon_matrix(ms[i], ms[j], strand='+', node_d=node_d)[0])

        mx = np.asarray([[1, 0, 0, 1, 0, 1, 0, 1, 1, 1]])
        # modified exon matrix
        self.assertFalse(compare_exon_matrix(ms[0], mx, strand='+', node_d=node_d)[0])

    def test_iterative_merge_transcripts(self):
        """No test yet"""
        pass # skip

    def test_get_fl_from_id(self):
        """Test get_fl_from_id(list_of_ids)"""
        ids = ["i2_HQ_sampleb92221|c242/f13p0/2019", "i2_HQ_sampleb92221|c24/f3p0/2019",
               "i2_HQ_sampleb92221|c103/f3p0/2371", "i2_HQ_sampleb92221|c107/f3p0/2368",
               "i2_HQ_sampleb92221|c1215/f2p0/2374", "i2_HQ_sampleb92221|c122/f4p0/2368"]
        expected_fl = 28
        self.assertEqual(get_fl_from_id(ids), expected_fl)

    def test_collapse_sam_records(self):
        """Test collapse_sam_records, which takes in a list of grouped sam records. and
        write collapsed gff records to good_gff_writer|bad_gff_writer. A collapsed
        gff record is 'good' if there are >= cov_threshold supportive sam records
        belonging to its group; otherwise, 'bad'.
        """
        test_name = "test_collapse_sam_records"
        good_gff_fn = op.join(_OUT_DIR_, test_name + ".good.gff.unfuzzy")
        bad_gff_fn = op.join(_OUT_DIR_, test_name + ".bad.gff.unfuzzy")
        group_fn = op.join(_OUT_DIR_, test_name + ".group.txt.unfuzzy")

        rmpath(good_gff_fn)
        rmpath(bad_gff_fn)
        rmpath(group_fn)

        records = _get_sam_groups()[0]["+"] # contains 10 sam records
        with CollapseGffWriter(good_gff_fn) as good_gff_writer, \
             CollapseGffWriter(bad_gff_fn) as  bad_gff_writer, \
             GroupWriter(group_fn) as group_writer:
            collapse_sam_records(records=records, cuff_index=0, cov_threshold=2,
                                 allow_extra_5exon=False, skip_5_exon_alt=True,
                                 good_gff_writer=good_gff_writer,
                                 bad_gff_writer=bad_gff_writer,
                                 group_writer=group_writer)

        def str_to_gffrecord(line):
            fields = line.strip().split('\t')
            print fields
            attributes = []
            for attr_tuple in fields[8].split(';'):
                if len(attr_tuple.strip()) == 0:
                    continue
                else:
                    fs = attr_tuple.strip().split(' ')
                    if len(fs) == 2:
                        attributes.append((fs[0], fs[1].replace('"', '')))
            return Gff3Record(seqid=fields[0], start=fields[3], end=fields[4],
                              type=fields[2], attributes=attributes)

        bad_gff_records = [str_to_gffrecord(line) for line in open(bad_gff_fn, 'r') if not line.startswith('##')]
        self.assertEqual(len(bad_gff_records), 0)

        good_gff_records = [str_to_gffrecord(line) for line in open(good_gff_fn, 'r') if not line.startswith('##')]

        self.assertEqual(len(good_gff_records), 4)
        self.assertEqual([(int(r.start), int(r.end), r.type, r.attributes['gene_id'], r.attributes['transcript_id']) for r in good_gff_records],
                         [(10711, 11641, 'transcript', "PB.0", "PB.0.1"),
                          (10711, 10791, 'exon', "PB.0", "PB.0.1"),
                          (10883, 11057, 'exon', "PB.0", "PB.0.1"),
                          (11435, 11641, 'exon', "PB.0", "PB.0.1"),
                         ])

    def test_Branch(self):
        """
        Test Branch and Branch.run.
        Note that fuzzy junctions are not merged.
        """
        test_name = "test_branch"
        good_gff_fn = op.join(_OUT_DIR_, test_name + ".good.gff.unfuzzy")
        bad_gff_fn = op.join(_OUT_DIR_, test_name + ".bad.gff.unfuzzy")
        group_fn = op.join(_OUT_DIR_, test_name + ".group.txt.unfuzzy")

        rmpath(good_gff_fn)
        rmpath(bad_gff_fn)
        rmpath(group_fn)

        b = Branch(isoform_filename=READS_DS, sam_filename=SORTED_GMAP_SAM,
                   cov_threshold=2, min_aln_coverage=0.99, min_aln_identity=0.95)

        b.run(allow_extra_5exon=True, skip_5_exon_alt=False,
              ignored_ids_fn=None,
              good_gff_fn=good_gff_fn,
              bad_gff_fn=bad_gff_fn,
              group_fn=group_fn)

        self.assertTrue(op.exists(good_gff_fn))
        self.assertTrue(op.exists(bad_gff_fn))
        self.assertTrue(op.exists(group_fn))

        std_good_gff_fn = op.join(SIV_STD_DIR, "test_branch", test_name + ".good.gff.unfuzzy")
        std_bad_gff_fn = op.join(SIV_STD_DIR, "test_branch", test_name + ".bad.gff.unfuzzy")
        std_group_fn = op.join(SIV_STD_DIR, "test_branch", test_name + ".group.txt.unfuzzy")

        print "Comparing %s and %s"  %  (good_gff_fn, std_good_gff_fn)
        self.assertTrue(filecmp.cmp(good_gff_fn, std_good_gff_fn))
        self.assertTrue(filecmp.cmp(bad_gff_fn, std_bad_gff_fn))
        self.assertTrue(filecmp.cmp(group_fn, std_group_fn))
