"""Test classes defined within pbtranscript.io.LA4IceReader."""
import unittest
import os.path as op
import hashlib
from pbtranscript.io import BLASRRecord, LA4IceReader
from test_setpath import DATA_DIR, OUT_DIR


t0_qaln="""ctatgagtaaat-atacta-gtata--a-atacga-ata-a--ata--tatctatagtatatataga--taga-tata-tag-agtatag-aata-aata-gta-ataatataagtata-tatatatagcta-a-tataatatactaata-gtaaaata-tactata-agtatatatactaagtatatatactatagtagtatataa-c-atagtagatata-ta-ta-atatatagtata-taatatatagtatagtata-taata-gtatatatactatagta-taataactatata-tgtatata-tata-gtatataggtataaata-gtataatatcccat-gtactct"""
t0_aaln="""|*|*||||***|*||||*|*|*|*|**|*||**||*|||*|**|||**|||||||||||||||*|||**||*|*||||*||**|*|||||*||||*||||*|||*|||||||||*||||*||||||||*|||*|*||||*|||||||*||*|||*|*||*|||||||*||||*|||||**||**|||||||*||**|||*|||||||*|*||||||*|||||*||*||*|||||||*||||*||*||||||*||*||||||*||*||*|||||||||||||||||*||*||*|||*|||*|*||||||*||||*|||*|||**||||*|||*|||||*|||***||*|||*|*|"""
t0_saln="""c-a-gagt-ggtgatacaacgcagagtacatgggatatatactatagatatctatagtatata-agaacta-agtataataatactatagtaatataatatgtatataatataactatagtatatata-ctatagtata-tatacta-taagtata-taatactatatagta-atata--aa--atatatagta--gta-tatataatctatagtatatatactaatatatatatactatagta-tatata-taaagtatacta-tatgtatatatactatagtaata-ta-cta-atagtatatatactataagta-ata--tatatataagtata-tatagtatagtaatat"""

t1_qaln="""aagcagtggt-at-caacgcagagtacatgggatatatactatagaataaatacatagtaatatataaatataatacgtatagtatataactaatagtaatataactatagtatata-tata-ctaataatagata-a--atag-ataaatagta-aa-tagatattactatg-agtaaatatactagtataaa-tacgaataaatatatctatagtatatatagatag-atatatagagtatagaataaatagtaataatataagtatatatatatag-ctaata-taatatactaatagtaaaatatactataagta-tatatac-taagtatatatactatagtagtatataacatag-tagatata-tataatatatagtata-taata-ta-tagta-ta-gtatata-ata-gtatatatactatagtataataactatatatgtatata-tatagtatat-aggtataaata-gtataatatccca-tg-ta"""
t1_aaln="""*||*||||||*||*||||||||||||||||||||||||||||||||*||**|*|*||*||*|||||||**|*|*||*|||||*||*|||*|||*|||||||||||*|||*|||||||*||||*|||*||*||*|||*|**||||*|||*|||*||*||*||*|||*||||||**|||||*||||**|*||||*|*||*|*|||*|||*||||||||||||||||**||**|||||||*|*|||||*|||*|||*|||*|*||||**|||||*||||||**|||*||*||||||||||||||||*|*|||||||||||||*||||||**||||||||||||*|||||||**|||||*|*|||*||*|||||*||||*||||||*||*|*||*||*||*|||||*||*|||||||*|||*|||||||||||||||||||*|||*||||*||**|||||*||*|*|||***|*|||****|**||*|*|||*||*|*||*|*"""
t1_saln="""cag-agtggtgatacaacgcagagtacatgggatatatactataga-ta--t-c-ta-tagtatataa-ga-actaagtataata-ata-cta-tagtaatataa-tat-gtatataatataacta-tagtatatatactatagtatatatactataagtatataatactatatagtaa-tata-aaatatatagta-gtatatata-atctatagtatatata-ctaatatatatatactatagtatatata-taa-agtata-ctatatgtatatatacta-tagtaatatactaatagtata-tatactataagtaatatatatataagtatatatagtatagta--atatatc-tagata-atatactatagtatata-ta-actagtagtaatagtaataagtatatatataagtatatatactatagtata-taa-tatagat-aatatagta-aatatcccatgta-ctctgcgtgtgata-cc-actgctt"""


class TEST_LA4ICEREADER(unittest.TestCase):
    """Test classes defined within pbtranscript.io.LA4IceReader."""
    def setUp(self):
        """Define las out file."""
        self.dataDir = DATA_DIR
        self.outDir = OUT_DIR
        self.las_out = op.join(self.dataDir, "test_LA4IceReader.las.out")
        # las_out was generated as follows:
        #   create reads.fasta containing four reads (t0_qaln, t0_saln, t1...)
        #   call 'ice_daligner.py reads.fasta reads.fasta out'
        #   copy the first two records in out/*.N1.las.out and stop signs to las_out
        self.t0 = BLASRRecord(
            qID=4,
            qLength=472, qStart=158, qEnd=472, qStrand='+',
            sID=3,
            sLength=461, sStart=0, sEnd=327, sStrand='+',
            score=-327, mapQV=None,
            qAln=t0_qaln, alnStr=t0_aaln, sAln=t0_saln,
            strand='+', identity=72.23)

        self.t1 = BLASRRecord(
            qID=4,
            qLength=472, qStart=0, qEnd=468, qStrand='+',
            sID=3,
            sLength=461, sStart=0, sEnd=461, sStrand='+',
            score=-461, mapQV=None,
            qAln=t1_qaln, alnStr=t1_aaln, sAln=t1_saln,
            strand='+', identity=75.03)


    def test_LA4IceReader(self):
        """Test LA4IceReader."""
        reader = LA4IceReader(self.las_out)
        reads = [x for x in reader]
        self.assertTrue(len(reads), 2)
        r0, r1 = reads

        self.assertTrue(r0 == self.t0)
        self.assertTrue(r1 == self.t1)

        f = op.join(self.outDir, "empty.las.out")
        with open(f, 'w') as writer:
            writer.write("+  +\n-  -\n")
        reads = [r for r in LA4IceReader(f)]
        self.assertTrue(len(reads) == 0)

