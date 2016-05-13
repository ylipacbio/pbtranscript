"""Test pClique, CliqueSubList."""
import unittest
import pbtranscript.ice.pClique as pClique
from networkx import Graph

class Test_pClique(unittest.TestCase):
    """Test pClique."""
    def test_maximal_cliques(self):
        """Test maximal_cliques."""
        G = Graph()
        G.add_edge('b','c')
        G.add_edge('b','d')
        G.add_edge('b','e')
        G.add_edge('b','f')
        G.add_edge('b','a')
        G.add_edge('a','c')
        G.add_edge('a','d')
        G.add_edge('a','e')
        G.add_edge('c','d')
        G.add_edge('c','e')
        G.add_edge('c','f')
        G.add_edge('c','g')
        G.add_edge('d','e')
        G.add_edge('d','g')
        G.add_edge('e','g')
        G.add_edge('f','g')
        # A clique of 'a, b, c, d, e' and some other edges.
        nodes = G.nodes()
        S, H = pClique.convert_graph_connectivity_to_sparse(G, nodes)
        i = nodes.index('a')
        tQ = pClique.grasp(S, H, 1, 5, i)
        c = [nodes[i] for i in tQ]
        print c
        self.assertTrue(set(c) == set(['a', 'b', 'c', 'd', 'e']))


if __name__ == "__main__":
    unittest.main()
