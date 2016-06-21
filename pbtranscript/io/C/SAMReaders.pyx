"""
Define classes to represent SAM records and read SAM records.
"""
from collections import namedtuple
import pysam
import logging

from pbtranscript.io.PbiBamIO import BamHeader

# __author__ = "etseng@pacificbiosciences.com"

__all__ = ["GMAPSAMReader", "GMAPSAMRecord", "iter_gmap_sam"]

Interval = namedtuple('Interval', ['start', 'end'])


def iter_cigar_string(cigar_string):
    """Iterate over a cigar string, each iter returns a tuple of
    (number, cigar_character) e.g., (10, 'M')"""
    num = cigar_string[0]
    for s in cigar_string[1:]:
        if str.isalpha(s):
            yield int(num), s
            num = ''
        else:
            num += s


def parse_cigar(cigar, offset):
    """
    Returns (segments, num_ins, num_del, num_mat_or_sub, query_start, query_end)
    where,
        segments is a list of Interval objects, of which start and end are
        genomic segment locations counting offset (e.g., offset = reference_start)
        num_ins, num_deletion: num of insertions, deletions w.r.t. reference
        num_mat_or_sub: num of matches or substitutions
        query_start, query_end: offset of query start and end position w.r.t.
        query read sequence.

    M - match
    I - insertion w.r.t. to ref
    D - deletion w.r.t. to ref
    N - skipped (which means splice junction)
    S - soft clipped
    H - hard clipped (not shown in SEQ)

    ex: 50M43N3D
    """
    segments = []
    query_start = 0
    cur_start = offset
    cur_end = offset
    first_thing = True
    q_aln_len = 0
    num_del = 0
    num_ins = 0
    num_mat_or_sub = 0
    for num, cigar_type in iter_cigar_string(cigar):
        if cigar_type == 'H' or cigar_type == 'S':
            if first_thing:
                query_start += num
        elif cigar_type == 'I':
            q_aln_len += num
            num_ins += num
        elif cigar_type == 'M':
            cur_end += num
            q_aln_len += num
            num_mat_or_sub += num
        elif cigar_type == 'D':
            cur_end += num
            num_del += num
        elif cigar_type == 'N': # junction, make a new segment
            segments.append(Interval(cur_start, cur_end))
            cur_start = cur_end + num
            cur_end = cur_start
        first_thing = False
    if cur_start != cur_end:
        segments.append(Interval(cur_start, cur_end))
    query_end = query_start + q_aln_len
    return segments, num_ins, num_del, num_mat_or_sub, query_start, query_end


class SAMflag(object):
    """SAMflag, can be constructed from int flag or by __init__."""
    __all__ = ('is_paired', 'strand', 'PE_read_num')
    def __init__(self, is_paired, strand, PE_read_num):
        self.is_paired = is_paired
        self.strand = strand
        self.PE_read_num = PE_read_num

    def __str__(self):
        return ("SAMflag(is_paired=%s, strand=%s, PE_read_num=%s)" %
                (self.is_paired, self.strand, self.PE_read_num))

    def __eq__(self, other):
        return (self.is_paired == other.is_paired and
                self.strand == other.strand and
                self.PE_read_num == other.PE_read_num)

    @classmethod
    def from_flag(cls, flag):
        """
        1 -- read is one of a pair
        2 -- alignment is one end of proper PE alignment          (IGNORE)
        4 -- read has no reported alignments                      (IGNORE)
        8 -- read is one of a pair and has no reported alignments (IGNORE)
        16 -- reverse ref strand
        32 -- other mate is aligned to ref strand
        64 -- first mate in pair
        128 -- second mate in pair
        256 -- not primary alignment

        Return: SAMflag
        """
        PE_read_num = 0
        strand = '+'
        if flag > 1024: #PCR or optical duplicate, should never see this...
            flag -= 1024
        if flag > 512: #not passing QC, should never see this
            flag -= 512
        if flag >= 256: #secondary alignment, OK to see this if option given in BowTie
            flag -= 256
        if flag >= 128:
            PE_read_num = 2
            flag -= 128
        if flag >= 64:
            PE_read_num = 1
            flag -= 64
        if flag >= 32:
            flag -= 32
        if flag >= 16:
            strand = '-'
            flag -= 16
        if flag >= 8:
            flag -= 8
        if flag >= 4:
            flag -= 4
        if flag >= 2:
            flag -= 2
        assert flag == 0 or flag == 1
        is_paired = (flag == 1)
        return SAMflag(is_paired, strand, PE_read_num)


class SAMRecordBase(object):
    """
    Base class for a SAMRecord.
    """
    def __init__(self):
        self.qID = None
        self.sID = None

        self.qStart = None
        self.qEnd = None

        self.sStart = None
        self.sEnd = None

        self.sLen = None
        self.qLen = None

        self.cigar = None
        self.flag = None

        self.segments = None

        self.num_nonmatches = None
        self.num_ins = None
        self.num_del = None
        self.num_mat_or_sub = None

        self.qCoverage = None
        self.sCoverage = None
        self.identity = None

    def __str__(self):
        msg =\
        """
        qID: {q}
        sID: {s}
        cigar: {c}
        sStart-sEnd: {ss}-{se}
        qStart-qEnd: {qs}-{qe}
        segments: {seg}
        flag: {f}

        coverage (of query): {qcov}
        coverage (of subject): {scov}
        alignment identity: {iden}
        """.format(q=self.qID, s=self.sID, seg=self.segments, c=self.cigar, f=self.flag,
                   ss=self.sStart, se=self.sEnd, qs=self.qStart, qe=self.qEnd,
                   iden=self.identity, qcov=self.qCoverage, scov=self.sCoverage)
        return msg

    def __eq__(self, other):
        def _cmp_(lhs, rhs):
            if lhs is None and rhs is None:
                return True
            else:
                return abs(lhs - rhs) <= 1e-6
        return (self.qID == other.qID and self.sID == other.sID and
                self.qStart == other.qStart and self.qEnd == self.qEnd and
                self.sStart == other.sStart and self.sEnd == other.sEnd and
                self.qLen == other.qLen and self.sLen == other.sLen and
                self.cigar == other.cigar and self.flag == other.flag and
                self.segments == other.segments and
                _cmp_(self.identity, other.identity) and
                _cmp_(self.qCoverage, other.qCoverage) and
                _cmp_(self.sCoverage, other.sCoverage))


class GMAPSAMRecord(SAMRecordBase):
    """Class represents GAMP SAM Record."""
    def __init__(self, alnseg, ref_len_dict=None, query_len_dict=None):
        """
        SAM files from GMAP have following optional fields:
            XS:A:+-   strand information
            NM:i:int  number of non-matches

        0. qID
        1. flag
        2. sID
        3. 1-based offset sStart
        4. mapping quality (ignore)
        5. cigar
        6. name of ref of mate alignment (ignore)
        7. 1-based offset sStart of mate (ignore)
        8. inferred fragment length (ignore)
        9. sequence (ignore)
        10. read qual (ignore)
        11. optional fields
        """
        super(GMAPSAMRecord, self).__init__()

        self.peer = alnseg

        self.qID = alnseg.query_name
        if alnseg.reference_id == -1: # means no match! STOP here
            self.sID = None
            return
        else:
            self.sID = alnseg.reference_name

        self.flag = SAMflag.from_flag(flag=alnseg.flag)
        self._flag_strand = self.flag.strand # serve as backup for debugging
        # Overwrite strand by XS:A:+- for GMAP SAM
        if alnseg.has_tag('XS') and alnseg.get_tag('XS') in ('+', '-'):
            self._flag_strand = self.flag.strand # serve as backup for debugging
            self.flag = SAMflag(is_paired=self.flag.is_paired,
                                strand=alnseg.get_tag('XS'),
                                PE_read_num=self.flag.PE_read_num)

        self.cigar = alnseg.cigarstring

        self.sStart = alnseg.reference_start
        (self.segments,
         self.num_ins, self.num_del, self.num_mat_or_sub,
         self.qStart, self.qEnd) = parse_cigar(cigar=self.cigar, offset=self.sStart)

        self.sEnd = self.segments[-1].end

        if ref_len_dict is not None:
            self.sCoverage = (self.sEnd - self.sStart) * 1. / ref_len_dict[self.sID]
            self.sLen = ref_len_dict[self.sID]

        if self.flag.strand == '-' and self.qLen is not None:
            self.qStart, self.qEnd = self.qLen - self.qEnd, self.qLen - self.qStart

        if self.qLen is not None:
            self.qCoverage = (self.qEnd - self.qStart) * 1. / self.qLen

        if query_len_dict is not None: # over write qLen and qCoverage, should be done LAST
            try:
                self.qLen = query_len_dict[self.qID]
            except KeyError: # HACK for blasr's extended qID
                k = self.qID.rfind('/')
                if k >= 0:
                    self.qLen = query_len_dict[self.qID[:self.qID.rfind('/')]]
                else:
                    self.qLen = query_len_dict[self.qID]
            self.qCoverage = (self.qEnd - self.qStart) * 1. / self.qLen

        if alnseg.has_tag('NM'):
            self.num_nonmatches = alnseg.get_tag('NM')
            self.identity = 1. - (self.num_nonmatches * 1. /
                                  (self.num_del + self.num_ins + self.num_mat_or_sub))

    @property
    def is_mapped(self):
        """Returns True if this SAM record is mapped (i.e., sID is neither None nor *)."""
        return self.sID is not None and self.sID != "*"


class GMAPSAMReader(object):
    """Class to read GMAP SAM file.,
    e.x.,
    with GMAPSAMReader("gmap-output.sam") as reader:
        for sam_record in reader:
            print sam_record
    """
    def __init__(self, filename, ref_len_dict=None, query_len_dict=None):
        self.filename = filename
        self._sam_file = pysam.Samfile(filename, 'r')
        self.header = BamHeader(headers=self._sam_file.header, ignore_pg=True)

        if ref_len_dict is not None:
            self.ref_len_dict = ref_len_dict
        else:
            self.ref_len_dict = self.header.referenceLengthsDict

        self.query_len_dict = query_len_dict

    def __iter__(self):
        return self

    def __str__(self):
        return "%s reading <%s>" % (self.__class__.__name__, self.filename)

    def next(self):
        """Return the next SAMRecord or stop."""
        return GMAPSAMRecord(alnseg=self._sam_file.next(),
                             ref_len_dict=self.ref_len_dict,
                             query_len_dict=self.query_len_dict)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._sam_file.close()


def iter_gmap_sam(sam_filename, query_len_dict,
                  min_aln_coverage, min_aln_identity, ignored_ids_writer):
    """
    Iterates over a SORTED GMAP SAM file, yields a collection of 
    GMAPSAMRecords that overlap by at least 1 base in the format:
    {'+': list of GMAPSAMRecords,
     '-': list of GMAPSAMRecords}

    Parameters:
      query_len_dict -- {query_read: query_read_len}, used for computing qCoverage
      ignored_ids_writer -- file handler writing ignored SAM record.
      min_aln_coverage -- ignore records of which qCoverage < min_aln_coverage 
      min_aln_identity -- ignore records of which identity < min_aln_identity
    """
    def sep_by_strand(records):
        """
        Returns {'+': list of + strand SAM records,
                 '-': list of - strand SAM record}
        """
        output = {'+':[], '-':[]}
        for r in records:
            output[r.flag.strand].append(r)
        return output

    def write_ignored_ids(msg):
        """Write ignored ids to ignored_ids_writer unless writer is None."""
        if ignored_ids_writer is not None:
            ignored_ids_writer.write(msg)

    records = None # holds the current set of records that overlap in coordinates

    with GMAPSAMReader(sam_filename, query_len_dict=query_len_dict) as reader:
        first_record = None
        for r in reader:
            if not r.is_mapped:
                write_ignored_ids("{0}\tUnmapped.\n".format(r.qID))
                #ignored_ids_writer.write("{0}\tUnmapped.\n".format(r.qID))
            elif r.qCoverage < min_aln_coverage:
                #ignored_ids_writer.write("{0}\tCoverage {1:.3f} too low.\n".format(r.qID, r.qCoverage))
                write_ignored_ids("{0}\tCoverage {1:.3f} too low.\n".format(r.qID, r.qCoverage))
            elif r.identity < min_aln_identity:
                #ignored_ids_writer.write("{0}\tIdentity {1:.3f} too low.\n".format(r.qID, r.identity))
                write_ignored_ids("{0}\tIdentity {1:.3f} too low.\n".format(r.qID, r.identity))
            else:
                first_record = r
                break
        if first_record is not None:
            records = [first_record]
        else:
            logging.warn("No valid records from %s!", sam_filename)
            return

        # The very first valid SAM record has been added to records,
        # continue to append remaining records and yield
        for r in reader:
            if r.sID == records[0].sID and r.sStart < records[-1].sStart:
                raise ValueError("SAM file %s is NOT sorted. ABORT!" % sam_filename)
            if r.qCoverage < min_aln_coverage:
                #ignored_ids_writer.write("{0}\tCoverage {1:.3f} too low.\n".format(r.qID, r.qCoverage))
                write_ignored_ids("{0}\tCoverage {1:.3f} too low.\n".format(r.qID, r.qCoverage))
            elif r.identity < min_aln_identity:
                #ignored_ids_writer.write("{0}\tIdentity {1:.3f} too low.\n".format(r.qID, r.identity))
                write_ignored_ids("{0}\tIdentity {1:.3f} too low.\n".format(r.qID, r.identity))
            elif r.sID != records[0].sID or r.sStart > max(x.sEnd for x in records):
                yield sep_by_strand(records)
                records = [r]
            else:
                records.append(r)
        yield sep_by_strand(records)
