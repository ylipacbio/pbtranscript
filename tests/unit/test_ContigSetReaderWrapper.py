"""Test classes defined within pbtranscript.io.ContigSetReaderWrapper."""
import unittest
import os.path as op
from pbtranscript.io import ContigSetReaderWrapper
from test_setpath import DATA_DIR, OUT_DIR


class TEST_ContigSetReaderWrapper(unittest.TestCase):
    """Test classes defined within pbtranscript.io.ContigSetReaderWrapper."""
    def setUp(self):
        """Define input and output file."""
        self.xmlfn = op.join(DATA_DIR, "test_ContigSetReaderWrapper.contigset.xml")
        self.fastafn = op.join(DATA_DIR, "reads_of_insert.fasta")
        self.fastqfn = op.join(DATA_DIR, "test_fastq_random_reader.fastq")

    def test_xml(self):
        """Test ContigSetReaderWrapper iterator over xml."""
        with ContigSetReaderWrapper(self.xmlfn) as reader:
            reads = [r for r in reader]
            self.assertEqual(len(reads), 25)

            self.assertEqual(reads[0].name, "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/43/ccs")
            self.assertTrue(reads[0].sequence[:].startswith("GTCCCAAATCCTGGGGAGTTCC"))

            self.assertEqual(reads[24].name, "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/54/ccs")
            self.assertTrue(reads[24].sequence[:].startswith("ACAAGTTTGACTTTGAAATCAGAG"))

        with ContigSetReaderWrapper(self.xmlfn) as reader:
            read = reader.next()
            self.assertEqual(read.name, "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/43/ccs")
            self.assertTrue(read.sequence[:].startswith("GTCCCAAATCCTGGGGAGTTCC"))

            out_fn = reader.consolidate(out_prefix=op.join(OUT_DIR, "test_ContigSetReaderWrapper_xml"))
            self.assertTrue(out_fn.endswith(".fasta"))
            self.assertTrue(len([r for r in ContigSetReaderWrapper(out_fn)]), 24)

    def test_fasta(self):
        """Test ContigSetReaderWrapper iterator over fasta."""
        with ContigSetReaderWrapper(self.fastafn) as reader:
            reads = [r for r in reader]
            self.assertEqual(len(reads), 22)

            self.assertEqual(reads[0].name, "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/43/ccs")
            self.assertTrue(reads[0].sequence[:].startswith("GTCCCAAATCCTGGGGAGTTCC"))

            self.assertEqual(reads[21].name, "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/970/ccs")
            self.assertTrue(reads[21].sequence[:].startswith("CTTACCAATGTGGGTCAGAT"))

        with ContigSetReaderWrapper(self.fastafn) as reader:
            read = reader.next()
            self.assertEqual(read.name, "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/43/ccs")
            self.assertTrue(read.sequence[:].startswith("GTCCCAAATCCTGGGGAGTTCC"))

            out_fn = reader.consolidate(out_prefix=op.join(OUT_DIR, "test_ContigSetReaderWrapper_fasta"))
            self.assertTrue(out_fn.endswith(".fasta"))
            self.assertTrue(len([r for r in ContigSetReaderWrapper(out_fn)]), 21)

    def test_fastq(self):
        """Test ContigSetReaderWrapper iterator over fastq."""
        with ContigSetReaderWrapper(self.fastqfn) as reader:
            reads = [r for r in reader]
            self.assertEqual(len(reads), 64)

            self.assertEqual(reads[0].name, "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/61232/4045_63_CCS")
            self.assertTrue(reads[0].sequence[:].startswith("ATTCTGAGGAGGGTCA"))

        with ContigSetReaderWrapper(self.fastqfn) as reader:
            read = reader.next()
            self.assertEqual(read.name, "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/61232/4045_63_CCS")
            self.assertTrue(read.sequence[:].startswith("ATTCTGAGGAGGGTCA"))

            out_fn = reader.consolidate(out_prefix=op.join(OUT_DIR, "test_ContigSetReaderWrapper_fastq"))
            self.assertTrue(out_fn.endswith(".fastq"))
            self.assertTrue(len([r for r in ContigSetReaderWrapper(out_fn)]), 63)

    def test_name_to_len_dict(self):
        expected = {"m131018_081703_42161_c100585152550000001823088404281404_s1_p0/61232/4045_63_CCS": 3982,
                    "m131018_081703_42161_c100585152550000001823088404281404_s1_p0/58475/3616_57_CCS":3559}
        d = ContigSetReaderWrapper.name_to_len_dict(self.fastqfn)
        for key in expected:
            self.assertEqual(expected[key], d[key])
