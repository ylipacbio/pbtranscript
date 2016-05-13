#!/usr/env python

"""
Define LA4IceReader which reads output of 'LA4Ice' as BLASRRecord.
"""

from pbtranscript.io.BLASRRecord import BLASRRecord

__author__ = 'etseng@pacificbiosciences.com'


class LA4IceReader(object):

    """
    Reader for reading alignments produced by
        'LA4Ice -m -i0 -w100000 -b0 -a {db} {las}'

    Example
        [r for r in LA4IceReader('*.las.out')][0] ==> this shows a BLASRRecord
    """

    def __init__(self, las_out_filename):
        self.file_name = las_out_filename
        self.f = self._open_file(las_out_filename)
        self._lineno = 0

    def _open_file(self, file_name):
        return open(file_name)

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        """Close *.las.out file."""
        self.f.close()

    def next(self):
        """
        Should be the printed out of running:
        LA4Ice -m -i0 -w100000 -b0 -a:{db} {las}
        """
        # first line is BLASR-like output
        # ex: 000000002 000002845 -1192 85.12 0 1855 3082 3082 0 2324 3516 3517 overlap
        # if is - strand, then strand=1, start=S, end=E means the sequence is
        # seq[S:E].reverse_complement()
        try:
            raw = self.f.readline().strip().split()
            self._lineno += 1
            if raw[0] == '+' and raw[1] == '+':  # FALCON-added EOF signature
                raise StopIteration
            qID = int(raw[0]) + 1  # convert to 1-based
            sID = int(raw[1]) + 1  # convert to 1-based
            score = int(raw[2])
            iden = float(raw[3])
            qStrand = int(raw[4])
            qStart = int(raw[5])  # 0-based
            qEnd = int(raw[6])
            qLen = int(raw[7])
            sStrand = int(raw[8])
            sStart = int(raw[9])
            sEnd = int(raw[10])
            sLen = int(raw[11])

            self.f.readline()  # blank line
            _qStart, qAln = self.f.readline().strip().split()
            self._lineno += 2
            # Liz: changed becuz new daligner has _qStart and _sStart both at 0-based
            assert ((qStrand == 0 and int(_qStart)  == qStart) or
                    (qStrand == 1 and int(_qStart)  == qEnd))
            alnStr = self.f.readline().strip()
            _sStart, sAln = self.f.readline().strip().split()[:2]
            self._lineno += 2
            # Liz: changed becuz new daligner has _qStart and _sStart both at 0-based
            assert ((sStrand == 0 and int(_sStart)  == sStart) or
                    (sStrand == 1 and int(_sStart)  == sEnd))
            return BLASRRecord(qID=qID, qLength=qLen,
                               qStart=qStart, qEnd=qEnd, qStrand=qStrand,
                               sID=sID, sLength=sLen,
                               sStart=sStart, sEnd=sEnd, sStrand=sStrand,
                               score=score, mapQV=None,
                               qAln=qAln, alnStr=alnStr, sAln=sAln,
                               identity=iden,
                               strand='+' if qStrand == sStrand else '-')
        except (IndexError, IOError, ValueError, AssertionError) as exc:
            raise ValueError("Unable to read %s line %d as LA4Ice output: %r." %
                             (self.file_name, self._lineno, exc))
