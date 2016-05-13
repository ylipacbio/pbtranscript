import unittest
import os.path as op
import filecmp
from cPickle import load
from pbtranscript.io import DazzIDHandler
from pbtranscript.Utils import mknewdir, execute
#from pbtranscript.io import *
from test_setpath import OUT_DIR, DATA_DIR, STD_DIR


class Test_DazzIDHandler(unittest.TestCase):
    """Test DazzIDHandler."""
    testName = "test_DazzIDHandler"
    def setUp(self):
        """Initialize."""
        self.inputDir  = op.join(DATA_DIR, self.testName)
        self.outDir    = op.join(OUT_DIR,  self.testName)
        self.stdoutDir = op.join(STD_DIR,  self.testName)
        self.fastaFileName = "test_DazzIDHandler.fasta"

        self.stdout_dazz_fasta = op.join(self.stdoutDir,
                                         self.fastaFileName[0:-6] + ".dazz.fasta")
        self.stdout_pickle = self.stdout_dazz_fasta + ".pickle"

        mknewdir(self.outDir)
        # Copy inputDir/test_DazzIDHandler.fasta to outDir.
        execute("cp %s %s" % (op.join(self.inputDir, self.fastaFileName),
                              op.join(self.outDir,   self.fastaFileName)))

    def test_convert_to_dazz_fasta(self):
        """convert input fasta to daligner-compatible fasta with ids."""
        fn = op.join(self.outDir, self.fastaFileName)
        prefix = fn[0:-6]

        expected_dazz_filename   = prefix + ".dazz.fasta"
        expected_pickle_filename = prefix + ".dazz.fasta.pickle"
        expected_db_filename = prefix + ".dazz.fasta.db"

        print "Testing DazzIDHandler(%s, converted=False)" % (fn)
        handler = DazzIDHandler(fn, converted=False)

        # check file names
        self.assertEqual(handler.dazz_filename,   expected_dazz_filename)
        self.assertEqual(handler.pickle_filename, expected_pickle_filename)
        self.assertEqual(handler.db_filename,     expected_db_filename)

        # converted files all exits
        self.assertTrue(op.exists(handler.dazz_filename))
        self.assertTrue(op.exists(handler.pickle_filename))

        # check converted daligner compatible fasta file
        self.assertTrue(filecmp.cmp(handler.dazz_filename, self.stdout_dazz_fasta))

        # check pickle
        self.assertTrue(filecmp.cmp(handler.pickle_filename, self.stdout_pickle))

        # check make_db
        print "Test DazzIDHandler.make_db"
        handler.make_db()
        self.assertTrue(op.exists(handler.db_filename))

        """Test self.read_dazz_pickle."""
        fn = op.join(self.stdoutDir, self.fastaFileName)
        print "Testing DazzIDHandler(%s, converted=True)" % fn
        handler = DazzIDHandler(fn, converted=True)
        handler.read_dazz_pickle()

        print "Loading from pickle %s" % (self.stdout_pickle)
        expected_dazz_mapping = load(open(self.stdout_pickle, 'r'))
        expected_keys = range(1, 13)

        self.assertTrue(expected_keys, handler.keys())
        for k in expected_keys:
            self.assertTrue(expected_dazz_mapping[k], handler.dazz_mapping[k])

        print "Testing DazzIDHandler.num_blocks"
        self.assertTrue(handler.num_blocks, 1)


class Test_DazzIDHandler_DataSet(unittest.TestCase):
    """Test DazzIDHandler while input is dataset."""
    testName = "test_DazzIDHandler_dataset"
    def setUp(self):
        """Initialize."""
        self.inputDir  = op.join(DATA_DIR, self.testName)
        self.outDir    = op.join(OUT_DIR,  self.testName)
        self.stdoutDir = op.join(STD_DIR,  self.testName)
        self.filename = "test_DazzIDHandler.contigset.xml"

        self.stdout_dazz_fasta = op.join(self.stdoutDir,
                                         self.filename[0:-4] + ".dazz.fasta")
        self.stdout_pickle = self.stdout_dazz_fasta + ".pickle"

        mknewdir(self.outDir)
        # Copy xml to outDir.
        #execute("cp %s %s" % (op.join(self.inputDir, self.filename),
        #                      op.join(self.outDir,   self.filename)))
        # Copy another.fasta to outDir
        #execute("cp %s %s" % (op.join(self.inputDir, "another.fasta"),
        #                      self.outDir))

    def test_convert_to_dazz_fasta(self):
        """convert input fasta to daligner-compatible fasta with ids."""
        prefix = op.join(self.outDir, self.filename)[0:-4]

        expected_dazz_filename   = prefix + ".dazz.fasta"
        expected_pickle_filename = prefix + ".dazz.fasta.pickle"
        expected_db_filename = prefix + ".dazz.fasta.db"

        print "Testing DazzIDHandler(xml, converted=False, dazz_dir)"
        handler = DazzIDHandler(op.join(self.inputDir, self.filename),
                                converted=False, dazz_dir=self.outDir)

        # check file names
        self.assertEqual(handler.dazz_filename,   expected_dazz_filename)
        self.assertEqual(handler.pickle_filename, expected_pickle_filename)
        self.assertEqual(handler.db_filename,     expected_db_filename)

        # converted files all exits
        self.assertTrue(op.exists(handler.dazz_filename))
        self.assertTrue(op.exists(handler.pickle_filename))

        # check converted daligner compatible fasta file
        self.assertTrue(filecmp.cmp(handler.dazz_filename, self.stdout_dazz_fasta))

        # check pickle
        self.assertTrue(filecmp.cmp(handler.pickle_filename, self.stdout_pickle))

        # check make_db
        print "Test DazzIDHandler.make_db"
        handler.make_db()
        self.assertTrue(op.exists(handler.db_filename))

        """Test self.read_dazz_pickle."""
        fn = op.join(self.stdoutDir, self.filename)
        print "Testing DazzIDHandler(%s, converted=True)" % fn
        handler = DazzIDHandler(fn, converted=True)
        handler.read_dazz_pickle()

        print "Loading from pickle %s" % (self.stdout_pickle)
        expected_dazz_mapping = load(open(self.stdout_pickle, 'r'))
        expected_keys = range(1, 16)

        self.assertTrue(expected_keys, handler.keys())
        for k in expected_keys:
            self.assertTrue(expected_dazz_mapping[k], handler.dazz_mapping[k])

        print "Testing DazzIDHandler.num_blocks"
        self.assertTrue(handler.num_blocks, 1)
