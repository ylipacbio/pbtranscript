# pylint: disable=R0921
# pylint: disable=W0223

"""
Defines classes for fetching data from a collection of original SMRTCell bam
files, and for writing bam files.
"""

__all__ = [ "PbiBamReader",
            "BamCollection",
            "BamHeader",
            "BamWriter"]


import os.path as op

import numpy as np
import pysam

from pbcore.util.Process import backticks
from pbcore.io import (DataSet, ReadSet, SubreadSet, ConsensusReadSet,
                       ContigSet, FastaReader, openDataFile, openDataSet)
from pbcore.io.dataset.DataSetMembers import Filters
from pbcore.io import BamAlignment

from pbtranscript.Utils import get_files_from_file_or_fofn, \
        guess_file_format, FILE_FORMATS


def _intersectRanges(r1, r2):
    b1, e1 = r1
    b2, e2 = r2
    b, e = max(b1, b2), min(e1, e2)
    return (b, e) if (b < e) else None


class BamZmw(object):

    """Represent zmws in bam file, like Zmw in pbcore.io"""

    def __init__(self, bamRecords, isCCS):
        if type(bamRecords) == list:
            assert(len(bamRecords) > 0)
            assert(isinstance(bamRecords[0], BamAlignment))
        else:
            assert(isinstance(bamRecords, BamAlignment))
            bamRecords = [bamRecords]
        assert(len(bamRecords) >= 1)
        self.bamRecords = bamRecords
        self._isCCS = isCCS
        if self._isCCS:
            assert(len(bamRecords) == 1)
        if not all([b.holeNumber == self.holeNumber for b in bamRecords]):
            raise ValueError("%s's input reads must come from same zmw." %
                              self.__class__.__name__)

    @property
    def holeNumber(self):
        """Return hole numer of this zmw."""
        return self.bamRecords[0].holeNumber

    @property
    def movieName(self):
        """Return movie name of this Zmw."""
        return self.bamRecords[0].qName.split("/")[0]

    @property
    def IsCCS(self):
        """True if this zmw has ccs reads; False otherwise."""
        return self._isCCS

    @property
    def zmwName(self):
        """Return zmw name in string."""
        return "/".join(self.bamRecords[0].qName.split("/")[0:2])

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__,
                             self.zmwName)

    @property
    def ccsRead(self):
        """Return ccs read of this zmw. None if this zmw only
        contains subreads."""
        if not self.IsCCS:
            return None
        assert(len(self.bamRecords) == 1)
        return BamCCSZmwRead(self, self.bamRecords[0], self.holeNumber)

    @property
    def subreads(self):
        """Return a list of subreads of this zmw."""
        if self.IsCCS:
            return []
        return [BamZmwRead(self, _record, self.holeNumber)
                for _record in self.bamRecords]

    def read(self, readStart=None, readEnd=None):
        """Return read of movie/zmw/readStart_readEnd.
           If either readStart or readEnd are None, return read
           of the zmw.
        """
        _record = None
        if readStart is None or readEnd is None:
            if len(self.bamRecords) != 1:
                # Reading full unrolled sequence in BamZmw is not supported
                raise ValueError("Must specify read start and end when" +
                                 "a BamZmw contains multiple bam records.")
            else:
                _record = self.bamRecords[0]
                readStart, readEnd = 0, _record.qLen
                if not self.IsCCS:
                    readStart, readEnd = _record.qStart, _record.qEnd
        else:
            _record = _findBamRecord(self.bamRecords, readStart, readEnd)
            if _record is None:
                raise IndexError("Invalid slice of zmw")
        return BamZmwRead(self, _record, self.holeNumber, readStart, readEnd)


def _findBamRecord(bamRecords, readStart, readEnd):
    """
    Given a list of BamAlignment object each of which represents a read,
    return the first record which covers interval [readStart, readEnd).
    """
    for bamRecord in bamRecords:
        if (bamRecord.qStart <= readStart and
            bamRecord.qEnd   >= readEnd):
            return bamRecord
    return None


class BamZmwRead(object):

    """Represent zmw read in bam files, like ZmwRead in pbcore.io"""
    __slots__ = ["zmw", "bamRecord", "holeNumber", "readStart", "readEnd"]

    def __init__(self, zmw, bamRecord, holeNumber,
                 readStart=None, readEnd=None):
        assert(type(zmw) == BamZmw)
        self.zmw = zmw
        self.bamRecord = bamRecord
        assert(type(self.bamRecord) == BamAlignment)
        assert(holeNumber == self.zmw.holeNumber)
        if readStart is not None:
            self.readStart = readStart
        else:
            self.readStart = self.qStart
        if readEnd is not None:
            self.readEnd = readEnd
        else:
            self.readEnd = self.qEnd
        if self.readStart < self.qStart or \
           self.readEnd > self.qEnd:
            raise ValueError("[%d, %d) not in %s" % (self.readStart,
                                                     self.readEnd,
                                                     self.bamRecord))

    @property
    def qStart(self):
        """
        Subread start pos in coordinate of zmw polymerase read, inclusive.
        """
        if self.bamRecord.isCCS:
            return 0
        else:
            return self.bamRecord.qStart

    @property
    def qEnd(self):
        """
        Subread end pos in coordinate of zmw polymerase read, exclusive.
        """
        if self.bamRecord.isCCS:
            return self.bamRecord.qLen
        else:
            return self.bamRecord.qEnd

    @property
    def movieName(self):
        """Return movie name of this read."""
        return self.zmw.movieName

    @property
    def holeNumber(self):
        """Return hole number of this read."""
        return self.zmw.holeNumber

    @property
    def readName(self):
        """Return zmw name string."""
        return "%s/%d_%d" % (self.zmw.zmwName,
                             self.readStart,
                             self.readEnd)

    @property
    def alignedSegment(self):
        """Return a pysam.calignmentfile.AlignedSegment
        object which is associated with this record.
        """
        return self.bamRecord.peer

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__,
                             self.readName)

    def __len__(self):
        return self.readEnd - self.readStart

    @property
    def _offsetBegin(self):
        """Read start offset."""
        return self.readStart - self.qStart

    @property
    def _offsetEnd(self):
        """Read end offset."""
        return self.readEnd - self.qStart

    def basecalls(self):
        """Return base calls of this zmw read in string."""
        return self.bamRecord.read(False)[self._offsetBegin: self._offsetEnd]

    def qv(self, _qvName):
        """Return QVs of this zmw read."""
        qvName = _qvName.lower()
        qvs = []
        try:
            if qvName == "deletionqv":
                qvs = self.bamRecord.DeletionQV(False)
            elif qvName == "insertionqv":
                qvs = self.bamRecord.InsertionQV(False)
            elif qvName == "mergeqv":
                qvs = self.bamRecord.MergeQV(False)
            elif qvName == "substitutionqv":
                qvs = self.bamRecord.SubstitutionQV(False)
            elif qvName == "substitutiontag":
                qvs = self.bamRecord.SubstitutionTag(False)
            elif qvName == "deletiontag":
                qvs = self.bamRecord.DeletionTag(False)
            elif qvName == "qualityvalue":
                if self.bamRecord.peer.qual is None:
                    raise AttributeError("No Quality Value")
                qvs = np.fromstring(self.bamRecord.peer.qual,
                                    dtype=np.uint8) - 33
            elif qvName == "ipd": # although not used by isoseq
                qvs = self.bamRecord.IPD(False)
            elif qvName == "pulsewidth": # although not used by isoseq
                qvs = self.bamRecord.PulseWidth(False)
            return qvs[self._offsetBegin: self._offsetEnd]
        except AttributeError:
            raise AttributeError("%s does not contain QV %s." %
                                 (self.readName, _qvName))

    def Clip(self, readStart, readEnd):
        """
        Clip this read to [readStart:readEnd), and return a copy of
        pysam.calignmentfile.AlignedSegment object.
        Assume that read.bamRecord.peer is an unmapped AlignedSegment.
        """
        new_query_name = "%s/%s/%d_%d" % (self.movieName, self.holeNumber,
                                          readStart, readEnd)
        if not (readStart >= self.readStart and readStart <= readEnd and
                readEnd <= self.readEnd):
            raise ValueError("Unable to clip subread %s from read %s." %
                             (new_query_name, self.readName))

        s, e = readStart - self.readStart, readEnd - self.readStart
        QV_TAGS = ["iq", "dq", "dt", "st", "sq", "mq", "ip", "pw"]

        # Create an unaligned pysam.AlignedSegment object.
        ret = pysam.AlignedSegment()
        ret.query_name = new_query_name

        peer = self.bamRecord.peer
        ret.query_sequence = peer.query_sequence[s:e]
        ret.flag = peer.flag

        assert(peer.reference_id == -1)
        assert(peer.reference_start == -1)
        assert(peer.cigartuples is None)
        assert(peer.next_reference_id == -1)
        assert(peer.next_reference_start == -1)
        assert(peer.template_length == 0)

        ret.reference_id = peer.reference_id
        ret.reference_start = peer.reference_start
        ret.cigar = []
        ret.next_reference_id = peer.next_reference_id
        ret.next_reference_start = peer.next_reference_start
        ret.template_length = peer.template_length

        ret.mapping_quality = peer.mapping_quality
        if peer.query_qualities is None:
            ret.query_qualities = None
        else:
            ret.query_qualities = pysam.array_to_qualitystring(peer.query_qualities)[s:e]

        tags = peer.tags[::]
        for index, (tag_name, tag_val) in enumerate(tags):
            if tag_name in QV_TAGS:
                if self.__len__() != len(tag_val):
                    raise ValueError ("%s's %s length %d ! = sequence length %d" %
                                      (peer.query_name, tag_name,
                                       len(tag_val), self.__len__()))
                tags[index] = (tag_name, tag_val[s:e])
            elif tag_name == 'qs':
                tags[index] = (tag_name, int(readStart))
            elif tag_name == 'qe':
                tags[index] = (tag_name, int(readEnd))

        ret.tags = tags
        return ret

    def __delitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__delitem__"))

    def __setitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__setitem__"))

    def __getitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__getitem__"))


class BamCCSZmwRead(BamZmwRead):

    """Represent Zmw ccs read."""
    @property
    def readName(self):
        return "%s/ccs" % self.zmw.zmwName


class ReaderAndRowNumber(object):

    """Associate BamReader objects with row number."""
    __slots__ = ["reader", "rowNumber", "bamRecord"]

    def __init__(self, reader, rowNumber):
        self.reader = reader
        self.rowNumber = rowNumber

    @property
    def bamRecord(self):
        """Return associated BamAlignment object."""
        return self.reader[self.rowNumber]


# FIXME this still assumes that individual .bam files contain reads from a
# single movie.  It works anyway if we specify a movie ID (as is done from
# BamCollection), but the code makes some assumptions about ZMW range.
class PbiBamReader(object):
    """
    Wrapper for providing keyed retrieval of reads.  Input dataset may be
    copied from another existing one, or read from a bam file or a fofn of bam
    files, assuming all reads come from the same SMRTCell, and are of the same
    type: either subreads or ccs.
    e.g., *.subreads.bam,
    *.1.subreads.bam, *.2.subreads.bam, *.3.subreads.bam, or
    *.1.ccs.bam, *.2.ccs.bam and *.3.ccs.bam.

        reader = PbiBamReader(bam)
        zmw8 = reader[8]
        subreads = zmw8.subreads
        ccsread  = zwm8.ccsRead
        read     = zwm8.read()
        subreads[0].readName
        subreads[0].basecalls()
        subreads[0].InsertionQV()
    """

    def __init__(self, dataset, *args, **kwds):
        if isinstance(dataset, DataSet):
            self._dataset = dataset
        else:
            self._dataset = openDataFile(*([dataset]+list(args)))
        if not isinstance(self._dataset, (SubreadSet, ConsensusReadSet)):
            raise TypeError("Input must be a SubreadSet or ConsensusReadSet.")
        self._movieName = kwds.get("movie_name", None)
        #if guess_file_format(self.bamFileNames) != FILE_FORMATS.BAM:
        #    raise IOError("%s can only read pacbio bam files or fofn" %
        #                  self.__class__.__name__)

        self._header = BamHeader(ignore_pg=True)
        self._indexedReaders = []
        self._hnranges = []

        for i_r, reader in enumerate(self._dataset.resourceReaders()):
            movie_ids = set([ rg.MovieName for rg in reader.readGroupTable ])
            if self._movieName is None:
                if len(movie_ids) > 1:
                    raise RuntimeError("Multiple movies found in %s" %
                        reader.filename)
                self._movieName = list(movie_ids)[0]
            else:
                for rg in reader.readGroupTable:
                    if rg.MovieName == self._movieName:
                        break
                else:
                    continue
            self._header.add(reader.peer.header)
            self._indexedReaders.append(i_r)
            # FIXME this seems risky
            _hnrange = (min(reader.holeNumber), max(reader.holeNumber))
            if any(_intersectRanges(_hnrange, hnr) is not None for hnr in self._hnranges):
                raise ValueError ("%s's zmw range [%d, %d] should not " %
                                  (reader.filename, _hnrange[0], _hnrange[1]) +
                                  "intersect with another file.")
            self._hnranges.append(_hnrange)

    def _get_reader(self, i_reader):
        """Return the i-th reader."""
        return self._dataset.resourceReaders()[i_reader]

    @property
    def movieName(self):
        """return movie name of smrtcells to read."""
        return self._movieName

    @property
    def readType(self):
        """CCS or SUBREAD?"""
        if isinstance(self._dataset, ConsensusReadSet):
            return "CCS"
        else:
            return "SUBREAD"

    def hn2reader(self, hn):
        """Return an IndexedBamReader in self._indexedReader which
        contains zmw hole number.
        self._indexedReaders should contain IndexedBamReaders for
        [*.1.subreads.bam, *.2.subreads.bam, *.3.subreads.bam]. Typical
        RS SMRTCell zmw ranges are:
        [(0, 54493), (54494, 108987), (108988, 170000).
        """
        for index, hnr in enumerate(self._hnranges):
            if hnr[0] <= hn and hnr[1] >= hn:
                return self._get_reader(self._indexedReaders[index])

    @property
    def header(self):
        """Merged bam header."""
        return self._header

    @property
    def IsCCS(self):
        """True if reads from a ccs bam."""
        return self.readType == "CCS"

    # FIXME note that this will probably blow away any attempted chunking by
    # ZMW
    def __iter__(self):
        """Iterate over ZMWs"""
        for _hn in sorted(set(self.holeNumber)):
            yield self[_hn]

    @property
    def holeNumber(self):
        """Return holeNumbers"""
        readers = [ self._get_reader(i) for i in self._indexedReaders ]
        return np.concatenate(tuple(ir.holeNumber for ir in readers))

    def subreads(self):
        """Iterate over all subreads of bam files."""
        # bam with mixed ccs reads and subreads not supported.
        if not self.IsCCS:
            for _hn in self.holeNumber:
                for subread in self._getitemScalar[_hn].subreads:
                    yield subread

    def ccsReads(self):
        """Iterate over all ccs reads of bam files."""
        if self.IsCCS:
            for _hn in self.holeNumber:
                yield self._getitemScalar(_hn).ccsRead

    def __len__(self):
        """Return total number of zmws."""
        return len(self.holeNumber)

    def _getitemScalar(self, hn):
        """Return BamZmw given hole number."""
        zmwName = "%s/%d" % (self.movieName, hn)
        hnreader = self.hn2reader(hn)

        if hn not in hnreader.holeNumber:
            raise ValueError("Could not access %s" % zmwName)
        else:
            return BamZmw(hnreader.readsByName(zmwName), isCCS=self.IsCCS)

    def __getitem__(self, holeNumbers):
        if (isinstance(holeNumbers, int) or
                issubclass(type(holeNumbers), np.integer)):
            return self._getitemScalar(holeNumbers)
        elif isinstance(holeNumbers, slice):
            raise NotImplementedError("%s.__getitem__(slice)" %
                                      self.__class__.__name__)
        elif (isinstance(holeNumbers, list) or
              isinstance(holeNumbers, np.ndarray)):
            if len(holeNumbers) == 0:
                return []
            else:
                entryType = type(holeNumbers[0])
                if entryType == int or issubclass(entryType, np.integer):
                    return [self._getitemScalar(r) for r in holeNumbers]
                elif entryType == bool or issubclass(entryType, np.bool_):
                    return [self._getitemScalar(r) for r in np.flatnonzero(holeNumbers)]
        raise TypeError("Invalid type for %s slicing" %
                        self.__class__.__name__)

    def __delitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__delitem__"))

    def __setitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__setitem__"))

    def close(self):
        """Close all bam readers."""
        self._dataset.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__, self.movieName)


class BamCollection(object):
    """
    Class representing a collection of `original` PacBio bam files.

    Can be initialized from a list of bam files, or an input.fofn
    file containing a list of bam files
    """

    def __init__(self, *args):
        if len(args) == 1:
            args = get_files_from_file_or_fofn(args[0])
        self._dataset = openDataFile(*args)
        # Implementation notes: find all the bam files, and group
        # them together by movieName
        movies = set()

        for bam in self._dataset.resourceReaders():
            for rg in bam.peer.header["RG"]: #readGroupTable:
                if rg['PL'] != "PACBIO":
                    raise IOError("Could not process file %s, " %
                        self.__class__.__name__ +
                        "input files must be orignal PacBio bam movies.")
            for rg in bam.readGroupTable:
                movies.add(rg.MovieName)
                assert rg.ReadType in ["CCS", "SUBREAD"]
        self._movieDataSets = {}
        for movie_id in movies:
            ds2 = self._dataset.copy()
            filt = Filters()
            filt.addRequirement(movie=[('==', movie_id)])
            ds2.addFilters(filt)
            self._movieDataSets[movie_id] = PbiBamReader(ds2,
                movie_name=movie_id)

    @property
    def movieNames(self):
        """Return movie names as a list of string."""
        return self._movieDataSets.keys()

    @property
    def readers(self):
        """Return all internal readers of this BamCollection."""
        return self._movieDataSets

    def _get_reader(self, movie_id):
        """Return reader for movie = movie_id."""
        ds = self._movieDataSets[movie_id]
        assert ds.movieName == movie_id
        return ds

    @property
    def header(self):
        """Return a BamHeader object which combines headers from
        all readers.
        """
        return BamHeader([reader.header for reader in self.readers.values()])

    def __getitem__(self, key):
        """
        Slice by movie name, zmw name, or zmw range name, using standard
        PacBio naming conventions.  Examples:

          - ["m110818_..._s1_p0"]             -> PbiBamReader
          - ["m110818_..._s1_p0/24480"]       -> Zmw
          - ["m110818_..._s1_p0/24480/20_67"] -> ZmwRead
          - ["m110818_..._s1_p0/24480/ccs"]   -> CCSZmwRead
        """
        if not isinstance(key, str):
            raise KeyError("key %s is not a string." % key)

        indices = key.rstrip("/").split("/")

        if len(indices) < 1:
            raise KeyError("Invalid slice of %s" % self.__class__.__name__)

        if len(indices) >= 1:
            result = self._get_reader(indices[0])
        if len(indices) >= 2:
            result = result[int(indices[1])]
        if len(indices) >= 3:
            if indices[2].lower() == "ccs":
                result = result.ccsRead
            else:
                start, end = map(int, indices[2].split("_"))
                result = result.read(start, end)
        return result

    #
    # Iterators over Zmw, ZmwRead objects
    #

    def __iter__(self):
        for reader in self.readers.values():
            for zmw in reader:
                yield zmw

    def reads(self):
        """Iterate over all reads"""
        for reader in self.readers.values():
            for read in reader.reads():
                yield read

    def subreads(self):
        """Iterate over all subreads."""
        for reader in self.readers.values():
            for read in reader.subreads():
                yield read

    def ccsReads(self):
        """Iterate over all ccs reads."""
        for reader in self.readers.values():
            for read in reader.ccsReads():
                yield read

    def __repr__(self):
        return "<%s of movies:\n%s>" % (self.__class__.__name__,
                "\n".join(self.movieNames))

    def __delitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__delitem__"))

    def __setitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__setitem__"))

    def __len__(self):
        """Return total number of zmws in all movies."""
        return sum([len(self.readers[m]) for m in self.movieNames])

    def close(self):
        """Close all readers."""
        for movie in self.movieNames:
            self.readers[movie].close()

    def __enter__(self):
        return self

    def __exit__(self):
        return self.close()


HDHEADERTAG = 'HD'
RGHEADERTAG = 'RG'
SQHEADERTAG = 'SQ'
PGHEADERTAG = 'PG'
HEADERTAGS  = [HDHEADERTAG, RGHEADERTAG, SQHEADERTAG, PGHEADERTAG]


class BamHeader(object):
    """
    Represent a bam header.
    """
    def __init__(self, headers=None, ignore_pg=False):
        """
        ignore_pg: whether or not to propogate @PG info.
        """
        self._dict = dict({tag:[] for tag in HEADERTAGS})
        self.ignore_pg = ignore_pg
        if headers is None:
            pass
        elif isinstance(headers, dict):
            assert(all([tag in headers for tag in HEADERTAGS]))
            self._dict = headers
        elif isinstance(headers, list):
            for _header in headers:
                self.add(_header)
        else:
            raise ValueError("<%s, __init__> does not support input type %s." %
                             (self.__class__.__name__, type(headers)))

    @property
    def headerLine(self):
        """Return the header line.
        e,g,. bh.headerLine = {'SO': 'UNKNOWN', 'VN': '1.5', 'pb': '3.0b5'}
        """
        return self._dict[HDHEADERTAG]

    @property
    def readGroups(self):
        """Return all read groups.
        e.g., bh.readGroups = [
        {'DS': 'READTYPE=CCS;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;BINDINGKIT=100236500;SEQUENCINGKIT=001558034;BASECALLERVERSION=2.3', 'ID': '9e743229','PL': 'PACBIO', 'PU': 'm131018_081703_42161_c100585152550000001823088404281404_s1_p0'},
        {'DS': 'READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;BINDINGKIT=100256000;SEQUENCINGKIT=100254800;BASECALLERVERSION=2.1.0.0.127824', 'ID': '0829d831','PL': 'PACBIO', 'PU': 'm140802_070303_42161_c110036822550000001823106706241502_s1_p0'}
        ]
        """
        return self._dict[RGHEADERTAG]

    @property
    def programGroups(self):
        """Return all program groups."""
        return self._dict[PGHEADERTAG]

    @property
    def referenceSequences(self):
        """Return all reference sequences."""
        return self._dict[SQHEADERTAG]

    @property
    def header(self):
        """Return this BamHeader as a dict"""
        return self._dict

    def containReadGroup(self, readGroupId):
        """Return whether or not contain this object contains
        a read group with ID=readGroupId.
        """
        return readGroupId in [rg['ID'] for rg in self.readGroups]

    def containReferenceSequence(self, refSeqId):
        """Return whether or not contain this object contains
        a reference sequence with ID=refSeqId.
        """
        return refSeqId in [refseq['SN'] for refseq in self.referenceSequences]

    def _addRG(self, rg):
        """Add a 'RG' entry to header, while rg['ID'] must be unique."""
        if not self.containReadGroup(rg['ID']):
            self.readGroups.append(rg)

    def _addPG(self, pg):
        """Add a 'PG' entry to header."""
        self.programGroups.append(pg)

    def _addSQ(self, sq):
        """Add a 'SQ' entry to header, while sq['SN'] must be unique."""
        if not self.containReferenceSequence(sq['SN']):
            self.referenceSequences.append(sq)

    def add(self, other):
        """Add/Merge another BamHeader or dict to this object."""
        _toadd = other
        if isinstance(other, BamHeader):
            _toadd = other.header
        else:
            if not isinstance(other, dict):
                raise TypeError("<%s, add> does not support type %s" %
                                (self.__class__.__name__, type(other)))

        if len(self.headerLine) == 0:
            self._dict[HDHEADERTAG] = _toadd[HDHEADERTAG]

        for rg in _toadd[RGHEADERTAG]:
            self._addRG(rg)

        if not self.ignore_pg:
            for pg in _toadd[PGHEADERTAG]:
                self._addPG(pg)

        for sq in _toadd[SQHEADERTAG]:
            self._addSQ(sq)

    def __repr__(self):
        maxn = 100
        return "<%s: %d RG, %d SN>\n" % (
                self.__class__.__name__,
                len(self.readGroups),
                len(self.referenceSequences)) + \
               "RGs [%s]\n" % (", ".join([rg['ID'] for rg in self.readGroups])[0:maxn]) + \
               "SNs [%s]\n" % (", ".join([sn['SN'] for sn in self.referenceSequences])[0:maxn])


class BamWriter(object):
    """
    Write unaligned bam records to a bam file.
    """

    def __init__(self, fname, header):
        self.filename = op.abspath(op.expanduser(fname))
        if isinstance(header, dict):
            self.peer = pysam.AlignmentFile(fname, "wb", header=header)
        elif isinstance(header, BamHeader):
            self.peer = pysam.AlignmentFile(fname, "wb", header=header.header)
        else:
            raise ValueError("<%s, __init__(fname, header)> " %
                             self.__class__.__name__ +
                             "header type must be either dict or BamHeader.")

    def close(self):
        """Close bam file."""
        self.peer.close()

    def write(self, record):
        """Write a record or a list of records.
        If record type is one of the followings,
            * pysam.calignmentfile.AlignedSegment
            * BamZmwRead
            * BamCCSZmwRead
        write this record.
        """
        # workaround for changes to the internal structure of pysam -
        # could probably be improved upon
        segment_class = getattr(pysam.calignmentfile, "AlignedSegment", None)
        if segment_class is None:
            segment_class = pysam.calignedsegment.AlignedSegment
        if isinstance(record, segment_class):
            self.peer.write(record)
        elif isinstance(record, BamCCSZmwRead) or \
             isinstance(record, BamZmwRead):
            self.peer.write(record.alignedSegment)
        else:
            raise ValueError("<%s, write()> does not support record type %s." %
                             (self.__class__.__name__, type(record)))

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__,
                             self.filename)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class CCSBamSequence(object):
    """
    Wrapper for a BAM record, used to substitute for a Fasta record.
    """
    def __init__(self, bam):
        self.name = bam.qname
        self.sequence = bam.seq


class CCSInput(object):
    """
    Wrapper class for handling multiple formats specifying CCS sequences.
    The old convention was to use .fasta, but we would like to be able to pass
    the classifier a ConsensusReadSet (i.e. .bam files) instead for use within
    pbsmrtpipe.
    """
    def __init__(self, file_name):
        self.file_name = file_name
        self._is_fasta = False
        self.ext = op.splitext(file_name)[1].upper()
        if self.ext in [".FA", ".FASTA"]:
            self._dataset = FastaReader(file_name)
            self._is_fasta = True
        elif self.ext == ".BAM":
            self._dataset = openDataFile(file_name)
        else: # either contigset.xml or consensusreadset.xml
            assert self.ext == ".XML"
            self._dataset = openDataSet(file_name)
            if isinstance(self._dataset, ContigSet):
                self._is_fasta = True

    def __iter__(self):
        for rec in self._dataset:
            if not self._is_fasta:
                rec = CCSBamSequence(rec.peer)
            yield rec

    def close(self):
        """Close all datasets."""
        self._dataset.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __len__(self):
        if not self._is_fasta:
            return len(self._dataset)
        else:
            if self.ext in [".FA", ".FASTA"]:
                return len([r for r in FastaReader(self.file_name)])
            else: # contigset
                n = 0
                for rr in self._dataset.resourceReaders():
                    n += len([r for r in rr])
                return n

    def __delitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__delitem__"))

    def __setitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__setitem__"))

    def __getitem__(self):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__getitem__"))

