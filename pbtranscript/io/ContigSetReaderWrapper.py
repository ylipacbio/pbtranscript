#!/usr/bin/env python

"""
class ContigSetReaderWrapper wraps ContigSet as reader
to read from FASTA, FASTQ, or xml files.

Note that: ContigSetReaderWrapper respects filters in
contigset xml files.

ContigSetReaderWrapper is provided to avoid the overhead
of using ContigSet when input is a FASTA or FASTQ file.
See difference between:
      [r for r in ContigSet("*.fasta")]
and
      [r for r in ContigSetWrapper("*.fasta")]
"""

import logging
import itertools
from pbcore.io import FastaReader, FastqReader, ContigSet

__author__ = 'yli@pacificbiosciences.com'

log = logging.getLogger(__name__)

__all__ = ["ContigSetReaderWrapper"]


class ContigSetReaderWrapper(object):

    """
    Wraps readers for input FASTA, FASTQ, ContigSet xml,

    <input>.fasta|fastq|xml accepted.

    e.g.:
        [r for r in ContigSetReaderWrapper(input_fn)]
    """
    FILE_TYPE = {"FA": "FASTA", "FASTA": "FASTA",
                 "FQ": "FASTQ", "FASTQ": "FASTQ",
                 "XML": "CONTIGSET"}
    def __init__(self, *input_filenames):
        """
        *input_filenames - input FASTA/FASTQ/ContigSet files
        """
        self.readers = self._open_files(*input_filenames)
        self.reader_index = 0
        self.it = self.readers[self.reader_index].__iter__()

    def get_file_type(self, input_filename):
        """Return file type: FASTA, FASTQ, CONTIGSET"""
        if not input_filename.rfind('.') >= 0:
            raise IOError("Could not recoginize file type of %s" % input_filename)
        else:
            suffix = input_filename[input_filename.rfind('.') + 1:].upper()
            return self.FILE_TYPE[suffix]

    def _open_files(self, *input_filenames):
        """Open file handers and return."""
        readers = []
        for fn in input_filenames:
            if self.get_file_type(fn) == "FASTA":
                readers.append(FastaReader(fn))
            elif self.get_file_type(fn) == "FASTQ":
                readers.append(FastqReader(fn))
            elif self.get_file_type(fn) == "CONTIGSET":
                readers.append(ContigSet(fn))
            else:
                raise IOError("Could not read %s as FASTA/FASTQ file." % fn)
        return readers

    def __iter__(self):
        iters = [reader.__iter__() for reader in self.readers]
        return itertools.chain(*iters)

    def __len__(self):
        errMsg = "%s.__len__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)

    def __delitem__(self, key):
        errMsg = "%s.__delitem__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)

    def __setitem__(self, key, value):
        errMsg = "%s.__setitem__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        """Close all readers."""
        for reader in self.readers:
            reader.close()

    def next(self):
        """Return the next FastaRecord or FastqRecord."""
        try:
            record = self.it.next()
            yield record
        except StopIteration:
            self.reader_index += 1
            if self.reader_index < len(self.readers):
                self.it = self.readers[self.reader_index].__iter__()
            else:
                raise StopIteration
