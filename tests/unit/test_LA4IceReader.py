"""Test classes defined within pbtranscript.io.LA4IceReader."""
import unittest
import os.path as op
import hashlib
from pbtranscript.io import BLASRRecord, LA4IceReader
from test_setpath import DATA_DIR, OUT_DIR


t0_qaln = """agaaacca-aga-agagaaagacaga-gaagagaacagaaaagaatacagaataggagagaaaaagagagagagacagagag-agaga-agaagag-ag-gaa-gaga-agata-gagagagaa-agaaa-gagga--ga-gagaagagagataa-acga-gagg-gagaagagagaagag-agaagagaacagaga-agagagaaaaaatgta-agagaagagagatagaaaa-aag-gagagat-a-gacgaa-g-aga-gagacagagaagag-agaacagatagaagagaatagagagaca-agaatgga-agcag-gaccgcg-gaacacgc--ggag"""
t0_aaln = """|||*||*|*|||*|*||||*||*|||*||||||||**||||*|||*|*||||*||*|||||||||||||||||||*|*|||**|||||*|||||||*||*|||*||||*|||*|*|||||||||*|||||*||*||**||*|||||*|||||*||*|*||*|||**||||*||||*|||||*|*||||*||*|||||*||*|||||**||*|*|*|*||||||*|||*|||*||*|||*||||*|**|*|||*||*|*|||*||||*||||||*||*||||*|||*|*|||||||*|||||||*|*||||*|*|*||*||*||*||*|*|||*||*|**||||"""
t0_saln = """agagac-agagagacagaacgaaagaagaagagaa--gaaa-gaagagagaagag-agagaaaaagagagagagagaaagaacagagacagaagaggagagaaagagagagaaaagagagagaagagaaacgaagacagaagagaacagaga-aagaagaagagaagaga-gagacaagagcaaaaga-aagagagagagcgagaaggaaag-agaaagaagaaagaaagacaagaagagagacaagaagacaaaaggagaagagagagagaa-aggagaaaaga-aaaagagaaaagagagagacagaa-gcagag-agagaacg-gagaagac-caaggag"""

t1_qaln = """cag-agtggtgatacaacgcagagtacatgggatatatactataga-ta--t-c-ta-tagtatataag-a-actaagtataata-ata-cta-tagtaatataa-tat-gtatataatataactatagta-tatatactatagtatatatactataagtatataatactatatagtaatataaaata-tatag-ta-gtatatata-atctatagtatatatac-taatatatatatactatagtatatata-taa-agtatac-tatatgtatatata-cta-tagtaatatactaatagtata-tatactataagtaatatatatataagtatatatagtatagtaatatat--c-tagata-atatactatagtatata-taactagtagtaatagtaataagtatatatataagtatatatactatagtata-taa-tatagata-atatagtaaa-tatcccatgtactctgc-gtgtgataccac-tg--ct-t"""
t1_aaln = """*||*||||||*||*||||||||||||||||||||||||||||||||*||**|*|*||*||*|||||||**|*|*||*|||||*||*|||*|||*|||||||||||*|||*|||||||*||||*|||*|*||*||*|||**||||*|||*|||*||*||*||*|||*||||||**|||||*|||*|*||*||||**||*|*|||*|||*||||||||||||||||**||**|||||||*|*|||||*|||*|||*|||*|*||||**|||||**|||||||*|||*||*||||||||||||||||*|*|||||||||||||*||||||**||||||||||||*|||||||*|||||**|*|||*||*|||||*||||*||||||*||**||*||*||*|||||*||*|||||||*|||*|||||||||||||||||||*|||*||||*||**|||||*||*|*|||***|*|||*|*****||*|*|||*|*|*||**||*|"""
t1_saln = """aagcagtggt-at-caacgcagagtacatgggatatatactatagaataaatacatagtaatatataaatataatacgtatagtatataactaatagtaatataactatagtatata-tata-cta-a-taatagata-aatag-ataaatagta-aa-tagatattactat-gagtaa-atatactagtataaatacgaataaatatatctatagtatatatagata-gatatatagagtatagaataaatagtaataatataagtatat--atatatagctaata-taatatactaatagtaaaatatactataagta-tatata-ctaagtatatatactatagtagtatataacatag-tagatata-tataatatatagtatataata-ta-tagta-ta-gtatata-ata-gtatatatactatagtataataactatatatgtatata-tatagtat-ataggta-taaatagtataatatcccatgtactct"""


class TEST_LA4ICEREADER(unittest.TestCase):
    """Test classes defined within pbtranscript.io.LA4IceReader."""
    def setUp(self):
        """Define las out file."""
        self.dataDir = DATA_DIR
        self.outDir = OUT_DIR
        self.las_out = op.join(self.dataDir, "test_LA4IceReader.las.out")
        self.t0 = BLASRRecord(
            qID=884,
            qLength=308, qStart=0, qEnd=308, qStrand='+',
            sID=472,
            sLength=410, sStart=74, sEnd=402, sStrand='+',
            score=-328, mapQV=None,
            qAln=t0_qaln, alnStr=t0_aaln, sAln=t0_saln,
            strand='+', identity=72.64)

        self.t1 = BLASRRecord(
            qID=1472,
            qLength=463, qStart=2, qEnd=463, qStrand='+',
            sID=480,
            sLength=489, sStart=0, sEnd=472, sStrand='+',
            score=-472, mapQV=None,
            qAln=t1_qaln, alnStr=t1_aaln, sAln=t1_saln,
            strand='+', identity=74.92)


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

