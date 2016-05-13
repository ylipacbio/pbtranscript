#!/usr/bin/env python

"""
class DazzIDHandler, which
(1) converts an arbitrary fasta file to a daligner-compatible fasta file.
(2) maintain mapping between read ids in daligner-compatible fasta file
    and its original name in input file.
(3) makes a dazz database so that input fasta file can run daligner later
"""

import logging
import os.path as op
from cPickle import load, dump
from pbcore.io import FastaWriter
from pbtranscript.io.ContigSetReaderWrapper import ContigSetReaderWrapper
from pbtranscript.Utils import execute, realpath, nfs_exists

__author__ = 'etseng@pacificbiosciences.com'

log = logging.getLogger(__name__)

class DazzIDHandler(object):

    """
    For any kind of input fasta, convert & maintain ID mapping

    <input>.fasta --> <input>.dazz.fasta, <input>.dazz.fasta.pickle

    If the input is fastq, automatically output as fasta since dalign only takes fasta

    By default, converted is False, so go ahead with converting.
    Otherwise, just read the mapping pickle.

    NOTE: Consistent with dazz's indexing, the IDs will be 1-based!!!!
    e.g.:
        cat <input>.fasta
        >movie/1001/200_300
        AA...GG
        >movie/2013/1_1000
        GC...GC

        cat <input>.dazz.fasta
        >prolog/1/0_100
        AA...GG
        >prolog/2/0_999
        GC...GC

        dump <input>.dazz.fasta.pickle
        {1: movie/1001/200_300, 2:movie/2013/1_1000}
    """
    dazz_movie_name = 'prolog'

    def __init__(self, input_filename, converted=False, dazz_dir=None):
        """
        input_filename - input FASTA/FASTQ/ContigSet file
        converted - whether or not input file has been converted to
                    daligner compatible FASTA file.
        dazz_dir - if None, save all dazz.fasta, dazz.pickle, db files
                  in the same directory as inputfile.
                  if a valid path, save all output files to dazz_dir.
        """
        self.dazz_dir = dazz_dir
        self.input_filename = realpath(input_filename)
        self.validate_file_type(self.input_filename)

        # index --> original sequence ID ex: 1 --> movie/zmw/start_end_CCS
        self.dazz_mapping = {}

        if converted and not nfs_exists(self.db_filename):
            log.warning(str(self.input_filename) +
                        " should have been converted to daligner-compatible" +
                        " format, but in fact it is not. Converting ...")
            converted = False

        if not converted:
            self.convert_to_dazz_fasta()
            self.make_db()
        else:
            self.read_dazz_pickle()

    def validate_file_type(self, filename):
        """Return input file type, either "fasta" or "fastq" or "xml",
        no other file type is allowed."""
        i = filename.rfind('.')
        file_suffix = filename[i + 1:].lower()
        if not file_suffix in ('fa', 'fasta', 'fastq', 'xml'):
            raise ValueError("%s " % (self.__class__.__name__) +
                             "can only take .fa, .fasta, .fq, .fastq" +
                             "or contigset.xml as input!")

    @property
    def dazz_filename(self):
        """Return converted daligner-compatible fasta file name.
        e.g., *.dazz.fasta
        """
        i = self.input_filename.rfind('.')
        filename = self.input_filename[:i] + '.dazz.fasta'
        if self.dazz_dir is not None:
            return op.join(self.dazz_dir, op.basename(filename))
        else:
            return filename

    @property
    def pickle_filename(self):
        """Return name of a pickle file which maps sequence names
        in the original fasta file to corresponding ids in
        daligner-compatibile fasta file.
        e.g., *.dazz.fasta.pickle
        """
        return self.dazz_filename + '.pickle'

    @property
    def db_filename(self):
        """Return file name of daligner database created taking
        daligner compatible fasta file as input.
        e.g., *.dazz.fasta.db
        """
        return self.dazz_filename + '.db'

    def convert_to_dazz_fasta(self):
        """
        Convert input fasta/fastq file to daligner-compatibe fasta with ids:
        <prefix>/<index>/0_<seqlen>

        Also write out mappings to pickle
        """
        log.debug("Converting %s to daligner compatible fasta %s.",
                  self.input_filename, self.dazz_filename)
        reader = ContigSetReaderWrapper(self.input_filename)

        with FastaWriter(self.dazz_filename) as f:
            i = 1
            for r in reader:
                f.writeRecord("{p}/{i}/0_{len}".format(p=self.dazz_movie_name,
                                                       i=i, len=len(r.sequence)),
                              r.sequence[:])
                self.dazz_mapping[i] = r.name
                i += 1

        reader.close()

        with open(self.pickle_filename, 'w') as f:
            dump(self.dazz_mapping, f)

    def read_dazz_pickle(self):
        """Read dazz mapping from pickle file."""
        log.debug("Reading daligner compatible fasta ids from pickle %s",
                  self.pickle_filename)
        try:
            with open(self.pickle_filename) as f:
                self.dazz_mapping = load(f)
        except IOError:
            log.debug("Could not open ", str(self.pickle_filename),
                      ", input fasta may not be converted.", )

    def make_db(self):
        """Make dazz database for input file.
        1. fasta2DB
        2. DBsplit
        3. get & store number of blocks
        *.dazz.fasta.db will be created.
        """
        log.debug("Making DAZZ database for %s.", self.dazz_filename)
        if not op.exists(self.dazz_filename):
            raise RuntimeError("%s hasn't been converted to daligner-compatible format." %
                               self.input_filename)
        if op.exists(self.db_filename):
            cmd = "DBrm %s" % self.dazz_filename
            execute(cmd=cmd)

        cmd = "fasta2DB %s %s " % (self.dazz_filename, self.dazz_filename)
        execute(cmd=cmd)

        cmd = "DBsplit -s200 %s" % self.dazz_filename
        execute(cmd)

    def keys(self):
        """Return all dazz ids in input FASTA/FASTQ."""
        return self.dazz_mapping.keys()

    @property
    def num_blocks(self):
        """Return number of blocks in DAZZ DB."""
        if not op.exists(self.db_filename):
            self.make_db()

        with open(self.db_filename) as f:
            f.readline()
            f.readline()
            x = f.readline().strip()
            if not x.startswith('blocks ='):
                raise ValueError("Could not get number of blocks in DAZZ DB %s"
                                 % self.db_filename)
        return int(x.split('=')[1])

    def __getitem__(self, key):
        """
        key should be a single integer from the id
        e.g.: prolog/1/0_1324 --> key should be 1
        """
        return self.dazz_mapping[key]
