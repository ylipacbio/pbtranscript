"""Test classes defined within pbtranscript.counting.CountUtils."""
import unittest
import os.path as op
from pbcore.util.Process import backticks
from pbtranscript.Utils import rmpath, mkdir
from pbtranscript.io import ReadStatReader, AbundanceReader
from pbtranscript.counting.CountingUtils import read_group_file, \
         output_read_count_FL, output_read_count_nFL, make_abundance_file
from test_setpath import DATA_DIR, OUT_DIR, SIV_DATA_DIR

_SIV_DIR_ = op.join(SIV_DATA_DIR, "test_counting")
_DAT_DIR_ = op.join(DATA_DIR, "test_counting")

GROUP_FN = op.join(_SIV_DIR_, "group.txt")

class TEST_CountUtils(unittest.TestCase):
    """Test functions of pbtranscript.counting.CountUtils."""
    def setUp(self):
        """Define input and output file."""
        pass

    def test_read_group_file(self):
        """Test read_group_file."""
        cid_info = read_group_file(group_filename=GROUP_FN,
                                   is_cid=True,
                                   sample_prefixes=None)

        self.assertEqual(len(cid_info[None].keys()), 846)
        self.assertEqual(cid_info[None]['i1_HQ_sampleb92221|c1030'], 'PB.5.6')
        self.assertEqual(cid_info[None]['i2_HQ_sampleb92221|c326'], 'PB.10.14')


        sample_prefixes = ['i0_HQ_sampleb92221', 'i1_HQ_sampleb92221', 'i2_HQ_sampleb92221']
        cid_info = read_group_file(group_filename=GROUP_FN,
                                   is_cid=True,
                                   sample_prefixes=sample_prefixes)
        self.assertEqual(len(cid_info[sample_prefixes[0]].keys()), 165)
        self.assertEqual(len(cid_info[sample_prefixes[1]].keys()), 297)
        self.assertEqual(len(cid_info[sample_prefixes[2]].keys()), 384)
        self.assertEqual(cid_info['i1_HQ_sampleb92221']['c1030'], 'PB.5.6')
        self.assertEqual(cid_info['i2_HQ_sampleb92221']['c326'], 'PB.10.14')

    def test_output_read_count_FL(self):
        """Test output_read_count_FL."""
        d = op.join(SIV_DATA_DIR, "test_make_abundance")
        bs = ["0to1kb_part0", "1to2kb_part0", "2to3kb_part0", "3to4kb_part0", "8to9kb_part0"]
        sample_prefixes = ["i%d_HQ_sample18ba5d" % i for i in (0, 1, 2, 3, 4)]
        pickles = [op.join(d, b, "cluster_out/output/final.pickle") for b in bs]
        group_filename = op.join(SIV_DATA_DIR, "test_make_abundance", "group.txt")
        prefix_pickle_tuples = zip(sample_prefixes, pickles)
        output_filename = op.join(OUT_DIR, "test_output_read_count_FL.read_stat.txt")
        output_mode = 'w'
        cid_info = read_group_file(group_filename=group_filename,
                                   is_cid=True, sample_prefixes=sample_prefixes)

        restricted_movies = None
        output_read_count_FL(cid_info=cid_info,
                             prefix_pickle_filename_tuples=prefix_pickle_tuples,
                             output_filename=output_filename,
                             output_mode=output_mode,
                             restricted_movies=restricted_movies)
        self.assertTrue(op.exists(output_filename))
        records = [r for r in ReadStatReader(output_filename)]
        self.assertEqual(len(records), 4712)

        expected_first = "m54006_160328_233933/39912051/31_505_CCS\t474\tY\tunique\tPB.1.1"
        expected_last = "m54006_160328_233933/47383436/629_57_CCS\t572\tY\tunmapped\tNA"

        self.assertEqual(str(records[0]), expected_first)
        self.assertEqual(str(records[-1]), expected_last)

        # Test with restricted movies
        output_filename = op.join(OUT_DIR, "test_output_read_count_FL.2.read_stat.txt")
        restricted_movies = ["m54006_160328_233933"]
        output_read_count_FL(cid_info=cid_info,
                             prefix_pickle_filename_tuples=prefix_pickle_tuples,
                             output_filename=output_filename,
                             output_mode=output_mode,
                             restricted_movies=restricted_movies)
        self.assertTrue(op.exists(output_filename))
        records = [r for r in ReadStatReader(output_filename)]
        self.assertEqual(len(records), 4712)
        self.assertEqual(str(records[0]), expected_first)
        self.assertEqual(str(records[-1]), expected_last)

    def test_output_read_count_nFL(self):
        """Test output_read_count_FL."""
        d = op.join(SIV_DATA_DIR, "test_make_abundance")
        bs = ["0to1kb_part0", "1to2kb_part0", "2to3kb_part0", "3to4kb_part0", "8to9kb_part0"]
        sample_prefixes = ["i%d_HQ_sample18ba5d" % i for i in (0, 1, 2, 3, 4)]
        pickles = [op.join(d, b, "cluster_out/output/map_noFL/nfl.all.partial_uc.pickle") for b in bs]
        group_filename = op.join(SIV_DATA_DIR, "test_make_abundance", "group.txt")
        prefix_pickle_tuples = zip(sample_prefixes, pickles)
        output_filename = op.join(OUT_DIR, "test_output_read_count_nFL.read_stat.txt")
        output_mode = 'w'
        cid_info = read_group_file(group_filename=group_filename,
                                   is_cid=True, sample_prefixes=sample_prefixes)

        restricted_movies = None
        output_read_count_nFL(cid_info=cid_info,
                             prefix_pickle_filename_tuples=prefix_pickle_tuples,
                             output_filename=output_filename,
                             output_mode=output_mode,
                             restricted_movies=restricted_movies)
        self.assertTrue(op.exists(output_filename))
        records = [r for r in ReadStatReader(output_filename)]
        self.assertEqual(len(records), 5703)

        expected_first = "m54006_160328_233933/11993579/0_2060_CCS\t2060\tN\tambiguous\tPB.5.4"
        expected_last = "m54006_160328_233933/23593293/0_1613_CCS\t1613\tN\tunmapped\tNA"
        self.assertEqual(str(records[0]), expected_first)
        self.assertEqual(str(records[-1]), expected_last)

        # Test with restricted movies
        output_filename = op.join(OUT_DIR, "test_output_read_count_nFL.2.read_stat.txt")
        restricted_movies = ["m54006_160328_233933"]
        output_read_count_nFL(cid_info=cid_info,
                             prefix_pickle_filename_tuples=prefix_pickle_tuples,
                             output_filename=output_filename,
                             output_mode=output_mode,
                             restricted_movies=restricted_movies)
        self.assertTrue(op.exists(output_filename))
        records = [r for r in ReadStatReader(output_filename)]
        self.assertEqual(len(records), 5703)

        expected_first = "m54006_160328_233933/11993579/0_2060_CCS\t2060\tN\tambiguous\tPB.5.4"
        expected_last = "m54006_160328_233933/37224924/0_2549_CCS\t2549\tN\tunmapped\tNA"
        self.assertEqual(str(records[0]), expected_first)
        self.assertEqual(str(records[-1]), expected_last)

    def test_make_abundance_file(self):
        """"""
        d = op.join(SIV_DATA_DIR, "test_make_abundance")
        read_stat_filename = op.join(d, "read_stat.txt")
        output_filename = op.join(OUT_DIR, "test_make_abundance_file.txt")
        make_abundance_file(read_stat_filename=read_stat_filename,
                            output_filename=output_filename,
                            given_total=None,
                            restricted_movies=None,
                            write_header_comments=True)
        print output_filename
        self.assertTrue(op.exists(output_filename))
        expected_first = "PB.1.1\t30\t30\t30.83\t6.3667e-03\t3.0367e-03\t3.0637e-03"
        expected_last = "PB.12.5\t16\t22\t23.20\t3.3956e-03\t2.2269e-03\t2.3052e-03"
        records = [r for r in AbundanceReader(output_filename)]
        self.assertEqual(len(records), 38)
        self.assertEqual(str(records[0]), expected_first)
        self.assertEqual(str(records[-1]), expected_last)
