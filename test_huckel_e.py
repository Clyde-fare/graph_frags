import unittest
from fragments import GraphFrags
from huckel import get_atoms_stability
from ase.io import read
import ASE_utils
from math import pi

__author__ = 'clyde'


class TestHuckel_e(unittest.TestCase):
    def setUp(self):

        test_mol = read('test1.xyz')
        test_mol.rotate('x', pi/2)
        test_mol.rotate('z', pi/2)

        m_test_mol = ASE_utils.to_molmod(test_mol)
        m_test_mol.set_default_graph()

        self.g = GraphFrags(m_test_mol)
        r = self.g.get_rings(10)[-1]
        self.g.set_base_frag([r])
        self.frag_it = self.g.gen_fragments()
        self.frags = list(self.frag_it)

        self.central_pair = {11, 12}

    def test_huckel_e(self):
        p_stab = get_atoms_stability(self.frags[-1], self.central_pair, self.g.master)
        self.assertAlmostEquals(p_stab, 25.518101628268539)


if __name__ == '__main__':
    unittest.main()