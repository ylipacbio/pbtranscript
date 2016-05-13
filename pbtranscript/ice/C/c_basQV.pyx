import sys
from pbcore.io import BasH5Reader
from pbtranscript.io import PbiBamReader
from cpython cimport bool
from libc.math cimport pow
from libcpp.deque cimport deque
from libcpp.vector cimport vector
cimport cython

ctypedef cython.int INTT

cdef inline double qv_to_prob(double qv):
    return pow(10, -qv / 10.)

cpdef precache_helper(char * bas_file, list seqids, list QV_names, dict qv_dict):
    cdef bool is_CCS
    cdef int s
    cdef int e
    cdef char strand
    cdef double qv

    # qv_dict = {} # seqid --> qv_name --> list of qv (transformed to prob)
    if (bas_file.endswith("h5")):
        bas = BasH5Reader(bas_file)
    elif bas_file.endswith("bam") or bas_file.endswith("xml"):
        bas = PbiBamReader(bas_file)
    else:
        raise IOError("Unable to precache QV for %s" % bas_file)

    # print >> sys.stderr, "worker precaching {0} ids from {1}".
    #format(len(seqids), bas_file)
    for seqid in seqids:  # must be in the same bas file
        seqid = seqid.split()[0]
        # in case there is extra information like fiveseen=1;threeseen=0;
        movie, hn, s_e = seqid.split('/')
        hn = int(hn)
        if s_e.endswith('_CCS'):
            is_CCS = True
            s, e = map(int, s_e.split('_')[:2])
        else:
            is_CCS = False
            s, e = map(int, s_e.split('_'))
        if s < e:
            strand = '+'
        else:
            s, e = e, s
            strand = '-'

        qv_dict[seqid] = {}
        for qv_name in QV_names:
            zmw = bas[hn]
            if is_CCS:
                if zmw.ccsRead is not None:
                    qvs = zmw.ccsRead.qv(qv_name)[s:e]
                else:  # this is for reads_of_insert w/ 0-passed
                    qvs = zmw.read(s, e).qv(qv_name)
            else:  # subread
                qvs = zmw.read(s, e).qv(qv_name)
            if strand == '-':
                qvs = qvs[::-1]
            qv_dict[seqid][qv_name] = [qv_to_prob(qv) for qv in qvs]
            del qvs
    del bas


def fastq_precache_helper(seqid, qvs, qv_dict):
    """
    similar to precache_helper except takes a QV quality 
    (ex: np array of [1, 13, 30])
    """
    qv_dict[seqid]['unsmoothed'] = [qv_to_prob(qv) for qv in qvs]


def maxval_per_window(list arr, int window_size):
    result = []
    maxval_per_window_helper(arr, window_size, result)
    return result


cdef maxval_per_window_helper(list arr, int window_size, list result):
    cdef deque[int] q
    cdef int i, j
    cdef double new_element
    cdef int len_arr
    cdef int w2

    len_arr = len(arr)
    w2 = window_size / 2

    # for 0th, compare arr[0] to arr[0+window_size/2]
    # for 1th, compare arr[0] to arr[1+window_size/2]
    # stop at arr[window_size/2]
    i = 0
    for j in xrange(1, w2 + 1):
        if arr[j] >= arr[i]:
            i = j
    q.push_back(i)
    result.append(arr[i])  # this is the max for 0-th element

    for i in xrange(-w2 + 1, len_arr - window_size + 1):
        # now looking at range [i, i+window_size)
        j = q.front()
        if j < i:
            q.pop_front()
        new_element = arr[i + window_size - 1]
        if q.empty():
            q.push_back(i + window_size - 1)
        elif new_element >= arr[j]:
            q.clear()
            q.push_back(i + window_size - 1)
        else:
            while not q.empty():
                j = q.back()
                q.pop_back()
                if arr[j] > new_element:
                    q.push_back(j)
                    break
            q.push_back(i + window_size - 1)
        j = q.front()
        result.append(arr[j])

    # finish the remainder
    for i in xrange(len_arr - w2, len_arr):
        j = q.front()
        while j < i - w2:
            q.pop_front()
            j = q.front()
        result.append(arr[j])

    q.clear()
