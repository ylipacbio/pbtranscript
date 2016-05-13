"""Test pbtranscript.Classifier."""

import unittest
import filecmp
import os.path as op

from pbcore.util.Process import backticks
from pbtranscript.tasks.TPickles import *
from pbtranscript.Utils import *
from test_setpath import DATA_DIR, OUT_DIR, STD_DIR, SIV_DATA_DIR

def chunk_task_i(cls, i):
    """Return an object of class cls with index i."""
    return cls(i, 'flnc_%s.fasta' % i, 'out_dir_%s' % i)


class Test_ChunkTask(unittest.TestCase):
    """Test ChunkTask."""
    def setUp(self):
        """Set up test data."""
        pass

    def test_write(self):
        """"""
        pass


class Test_ChunkTasksPickle(unittest.TestCase):
    """Test ChunkTasksPickle."""
    def setUp(self):
        """Set up test data"""
        self.NUM = 4
        self.out_dir = op.join(OUT_DIR, "test_ChunkTasksPickle")
        mkdir(self.out_dir)

        self.chunk_tasks = [chunk_task_i(cls=ChunkTask, i=i) for i in range(0, self.NUM)]
        self.append_chunk_task = chunk_task_i(cls=ChunkTask, i=self.NUM)
        self.out_pickle_fns = [op.join(self.out_dir, '%s.pickle' % i)
                               for i in range(0, self.NUM + 1)]

    def test_read_write(self):
        """Test write."""
        fn = op.join(self.out_dir, "test_write.pickle")
        p = ChunkTasksPickle(self.chunk_tasks)
        p.write(fn)
        q = ChunkTasksPickle.read(fn)
        self.assertEqual(len(q), self.NUM)

        self.assertTrue(all([isinstance(r, ChunkTask) for r in p]))

        p.append(self.append_chunk_task)
        self.assertEqual(len(p), self.NUM + 1)

        p.spawn_pickles(self.out_pickle_fns)
        tmp = self.chunk_tasks + [self.append_chunk_task]
        for i, fn in enumerate(self.out_pickle_fns):
            self.assertTrue(op.exists(fn))
            self.assertEqual(len(ChunkTasksPickle.read(fn)), 1)
            self.assertEqual(ChunkTasksPickle.read(fn)[0], tmp[i])

    def test_sorted_by_attr(self):
        """Test sorted by attribute"""
        unsorted_indices = (10, 1, 5, 4, 10, 4, 5, 2)
        chunk_tasks = [chunk_task_i(cls=ChunkTask, i=i) for i in unsorted_indices]
        sorted_chunk_tasks = [chunk_task_i(cls=ChunkTask, i=i) for i in sorted(unsorted_indices)]
        p = ChunkTasksPickle(chunk_tasks)
        p.sorted_by_attr(attr='cluster_bin_index')
        self.assertEqual(p.chunk_tasks, sorted_chunk_tasks)

    def test_sort_and_group_tasks(self):
        """Test sort_and_group_tasks"""
        d = op.join(SIV_DATA_DIR, "test_tool_contract_chunks")
        p_fn = op.join(d, "cluster_chunks.pickle")
        print p_fn
        p = ChunkTasksPickle.read(p_fn)
        groups = p.sort_and_group_tasks(max_nchunks=10)
        print 'groups=%s' % groups
        expected_groups = [[0], [1], [2], [3], [4]]
        self.assertEqual(groups, expected_groups)

        groups = p.sort_and_group_tasks(max_nchunks=1)
        print 'groups=%s' % groups
        expected_groups = [[0,1,2,3,4]]
        self.assertEqual(groups, expected_groups)

        p_fn = op.join(d, "partial_chunks.pickle")
        print p_fn
        p = ChunkTasksPickle.read(p_fn)
        groups = p.sort_and_group_tasks(max_nchunks=8)
        print 'groups=%s' % groups
        expected_groups = [[0,8], [1,9], [2], [3], [4], [5], [6], [7]]
        self.assertEqual(groups, expected_groups)
        groups = p.sort_and_group_tasks(max_nchunks=1)
        print 'groups=%s' % groups
        expected_groups = [[0,1,2,3,4,5,6,7,8,9]]
        self.assertEqual(groups, expected_groups)

        p_fn = op.join(d, "polish_chunks.pickle")
        print p_fn
        p = ChunkTasksPickle.read(p_fn)
        groups = p.sort_and_group_tasks(max_nchunks=10)
        print 'groups=%s' % groups
        expected_groups = [[0, 10], [1,11], [2,12], [3], [4], [5], [6], [7], [8], [9]]
        self.assertEqual(groups, expected_groups)
        groups = p.sort_and_group_tasks(max_nchunks=1)
        print 'groups=%s' % groups
        expected_groups = [[0,1,2,3,4,5,6,7,8,9,10,11,12]]
        self.assertEqual(groups, expected_groups)
