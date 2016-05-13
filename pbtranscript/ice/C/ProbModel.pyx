"""Define Probability Models, including ProbFromQV and ProbFromModel."""
# distutils: language = c++
# distutils: sources = ProbModel.cpp
from pbtranscript.io.BasQV import basQVcacher, fastqQVcacher
from pbtranscript.io.ContigSetReaderWrapper import ContigSetReaderWrapper
from pbcore.io import FastqReader
from libc.math cimport log

DEFAULT_WINDOW_SIZE = 3 # default window size changed from 5 to 3.

class ProbFromFastq:
    """
    Probability model constructed from Fastq files using 
    a single QV for everything
    """
    def __init__(self, fastq_filename,
                 prob_threshold=.1, window_size=DEFAULT_WINDOW_SIZE):
        self.qver = fastqQVcacher()
        self.fastq_filename = fastq_filename
        self.seqids = []
        self.prob_threshold = prob_threshold
        self.window_size = window_size
        self.full_prob = None

        self.add_seqs_from_fastq(fastq_filename)

    def get_smoothed(self, qID, qvname, position=None):
        """
        Get smoothed QV of read=qID, type=qvname, position=position.
        In reality, qvname is ignored.
        """
        return self.qver.get_smoothed(qID, qvname, position)

    def get(self, qID, qvname, position=None):
        """
        Get QV of read=qID, type=qvname, position=position.
        In reality, qvname is ignored.
        """
        return self.qver.get(qID, qvname, position)
    
    def add_seqs_from_fastq(self, fastq_filename, smooth=True):
        """Add sequence ids from a fastq file."""
        self.qver.precache_fastq(fastq_filename)
        newids = [r.name.split()[0] for r in FastqReader(fastq_filename)]
        self.seqids += newids
        if smooth:
            self.qver.presmooth(newids, self.window_size, fastq_filename)

    def get_mean(self, qID, qvname):
        """Return mean QV of read=qID, type=qvname"""
        return self.qver.get_mean(qID, qvname)

    def remove_ids(self, ids):
        """Remove ids from self.seqids."""
        for _id in ids:
            self.seqids.remove(_id)
            del self.qver.qv[_id]

    def calc_prob_from_aln(self, qID, qStart, qEnd, fakecigar):
        """
        aln_record is a pysam.AlignedRead
        """
        prob_err = self.qver.get(qID, None)
        return calc_aln_log_prob2(prob_err, len(prob_err), list(fakecigar), 
                                  qStart, qEnd)
        

class ProbFromQV:

    """
    Probability model constructed from FOFN files using
    quality values.
    """

    def __init__(self, input_fofn, fasta_filename=None,
                 prob_threshold=.1, window_size=DEFAULT_WINDOW_SIZE):
        self.qver = basQVcacher()
        self.input_fofn = input_fofn
        self.seqids = []
        self.prob_threshold = prob_threshold
        self.window_size = window_size

        if self.input_fofn.endswith(".consensusreadset.xml"):
            self.qver.add_bash5(self.input_fofn)
        else:
            with open(self.input_fofn) as f:
                for line in f:
                    self.qver.add_bash5(line.strip())

        if fasta_filename is not None:
            self.add_seqs_from_fasta(fasta_filename)

    def get_smoothed(self, qID, qvname, position=None):
        """
        Get smoothed QV of read=qID, type=qvname, position=position.
        """
        return self.qver.get_smoothed(qID, qvname, position)

    def get(self, qID, qvname, position=None):
        """
        Get QV of read=qID, type=qvname, position=position.
        """
        return self.qver.get(qID, qvname, position)

    def get_mean(self, qID, qvname):
        """Return mean QV of read=qID, type=qvname."""
        return self.qver.get_mean(qID, qvname)

    def add_seqs_from_fasta(self, fasta_filename, smooth=True):
        """Add sequence ids from a fasta file."""
        with ContigSetReaderWrapper(fasta_filename) as reader:
            newids = [r.name.split()[0] for r in reader]
        self.add_ids_from_fasta(newids)

    def add_ids_from_fasta(self, newids):
        """Add sequence ids."""
        self.qver.precache(newids)
        self.seqids += newids
        self.qver.presmooth(newids, self.window_size)

    def remove_ids(self, ids):
        """Remove ids from self.seqids."""
        for _id in ids:
            self.seqids.remove(_id)
            del self.qver.qv[_id]

    def calc_prob_from_aln(self, qID, qStart, qEnd, fakecigar):
        """
        Calculate substitution/insertion/deletion probabilities.
        """
        prob_sub = self.qver.get(qID, 'SubstitutionQV')
        prob_ins = self.qver.get(qID, 'InsertionQV')
        prob_del = self.qver.get(qID, 'DeletionQV')
        return calc_aln_log_prob(prob_sub, prob_ins, prob_del,
                                 len(prob_del), list(fakecigar), 
                                 qStart, qEnd)


cdef double calc_aln_log_prob(list prob_sub, list prob_ins,
                              list prob_del, int n,
                              list fakecigar, int qStart, int qEnd):
    """
    Calculate log probabilities from alignment cigar strings using Cython.
    """
    cdef int i, cur_q_pos
    cdef double score, tmp, one_three

    one_three = log(1 / 3.)
    cur_q_pos = qStart
    score = 0.
    for x in fakecigar:
        if x == 'M':
            tmp = 1 - prob_sub[cur_q_pos] - \
                prob_ins[cur_q_pos] - prob_del[cur_q_pos]
            # sanity check...sometimes this can have < 0 prob,
            # so assign it a small prob like 0.001
            if tmp <= 0:
                tmp = 0.001
            score += log(tmp)
            cur_q_pos += 1
        elif x == 'S':
            score += log(prob_sub[cur_q_pos]) + one_three
            cur_q_pos += 1
        elif x == 'I':
            score += log(prob_ins[cur_q_pos]) + one_three
            cur_q_pos += 1
        else:  # x == 'D', don't advance qpos
            score += log(prob_del[cur_q_pos])
    assert cur_q_pos == qEnd
    return score


cdef double calc_aln_log_prob2(list prob_err, int n,
                               list fakecigar, int qStart, int qEnd):
    """
    Calculate log probabilities from alignment cigar strings using Cython.
    """
    cdef int i, cur_q_pos
    cdef double score, tmp, one_three

    cdef int max_pos = len(prob_err)
    one_three = log(1 / 3.)
    cur_q_pos = qStart
    score = 0.
    for x in fakecigar:
        if cur_q_pos >= max_pos: 
            # ToDo: this is a bug caused by daligner coordinates issues,
            # FIX LATER instead of just returning as done here
            return score
        if x == 'M':
            tmp = 1 - prob_err[cur_q_pos] 
            # sanity check...sometimes this can have < 0 prob,
            # so assign it a small prob like 0.001
            if tmp <= 0:
                tmp = 0.001
            score += log(tmp)
            cur_q_pos += 1
        elif x == 'S':
            score += log(prob_err[cur_q_pos]) + one_three
            cur_q_pos += 1
        elif x == 'I':
            score += log(prob_err[cur_q_pos]) + one_three
            cur_q_pos += 1
        else:  # x == 'D', don't advance qpos
            score += log(prob_err[cur_q_pos])
    assert cur_q_pos == qEnd
    return score


class fakeQVer(object):

    """
    Used by ProbFromModel to support the fake .get and .getsmoothed
    """

    def __init__(self, r_mis, r_ins, r_del):
        self.r_mis = r_mis
        self.r_ins = r_ins
        self.r_del = r_del
        self.err_mean = r_mis + r_ins + r_del
        self.d = {'DeletionQV': r_del, 'InsertionQV': r_ins, 'SubstitutionQV': r_mis}

    def get(self, qID, qvname, position=None):
        """Fake get()."""
        return self.d[qvname]

    def get_smoothed(self, qID, qvname, position=None):
        """Fake get_smoothed."""
        return self.d[qvname]

    def get_mean(self, qID, qvname):
        """Return fake mean QV."""
        return self.err_mean


class ProbFromModel:

    """Probability model from fixed indel/substitution rates."""

    def __init__(self, r_mis, r_ins, r_del):
        self.r_mis = r_mis
        self.r_ins = r_ins
        self.r_del = r_del
        self.r_mat = 1 - r_mis - r_ins - r_del
        assert self.r_mat > 0
        self.qver = fakeQVer(r_mis, r_ins, r_del)

    def add_seqs_from_fasta(self, fasta_filename, smooth=True):
        """dummy, nothing to do"""
        pass

    def remove_ids(self, ids):
        """dummy, nothing to do"""
        pass

    def get_qver_smoothed(self, *args):
        """Fake get_qver_smoothed."""
        return 0.

    def get(self, qID, qvname, position=None):
        """Return QV for qID, qvname, position."""
        return self.qver.get(qID, qvname, position)

    def get_smoothed(self, qID, qvname, position=None):
        """Get smoothed QV for qID, qvname, position."""
        return self.qver.get_smoothed(qID, qvname, position=position)

    def get_mean(self, seqid, qv_name):
        """Return mean QV for qID, qvname."""
        return self.qver.get_mean(seqid, qv_name)

    def calc_prob_from_aln(self, qID, qStart, qEnd, fakecigar):
        """Calculate probability from an alingment."""
        prob_mat = log(self.r_mat)
        prob_sub = log(self.r_mis)
        prob_ins = log(self.r_ins)
        prob_del = log(self.r_del)

        score = 0
        for x in fakecigar:
            if x == 'M':
                score += prob_mat
            elif x == 'S':
                score += prob_sub
            elif x == 'I':
                score += prob_ins
            else:  # x == 'D', don't advance qpos
                score += prob_del
        return score
