import unittest
from fragments import GraphFrags
from utils import get_unique
from ase.io import read
import ASE_utils
from math import pi

__author__ = 'clyde'


class TestGraphFrags(unittest.TestCase):
    def setUp(self):
        test_mol = read('test1.xyz')
        test_mol.rotate('x', pi/2)
        test_mol.rotate('z', pi/2)
        m_test_mol = ASE_utils.to_molmod(test_mol)
        m_test_mol.set_default_graph()
        self.mol_graph = GraphFrags(m_test_mol)

    def test_0(self):
        """Tests first shell of neighboring rings correctly identified"""

        initial_ring = self.mol_graph.get_rings(10)[-1]
        self.mol_graph.set_base_frag([initial_ring])
        first_shell = self.mol_graph.get_allowed_rings(1)

        assert first_shell == (frozenset({0, 1, 8, 9, 10, 11}),
                               frozenset({9, 10, 17, 18, 19, 20}),
                               frozenset({1, 2, 3, 4, 11, 12}),
                               frozenset({4, 5, 12, 13, 14, 15}),
                               frozenset({19, 20, 21, 22, 25, 26}),
                               frozenset({13, 14, 21, 22, 23, 24}))

    def test_1a(self):
        """Tests the number of fragments found matches a known case"""

        initial_ring = self.mol_graph.get_rings(10)[-1]
        self.mol_graph.set_base_frag([initial_ring])
        frag_it = self.mol_graph.gen_fragments()

        assert len(list(frag_it)) == 112


    def test_1b(self):
        """Tests no redundant fragments are produced"""

        initial_ring = self.mol_graph.get_rings(10)[-1]
        self.mol_graph.set_base_frag([initial_ring])
        frag_it = self.mol_graph.gen_fragments()

        list_frags = list(frag_it)
        assert len(get_unique(list_frags)) == len(list_frags)

if __name__ == '__main__':
    unittest.main()