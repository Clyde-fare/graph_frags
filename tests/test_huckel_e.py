import unittest
import ASE_utils
import huckel
from fragments import GraphFrags
from ase.io import read
from ase.structure import molecule
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
        p_stab = huckel.get_atoms_stability(self.frags[-1], self.central_pair, self.g.master)
        self.assertAlmostEquals(p_stab, 25.518101628268539)

    def test_huckel_e2(self):
        ethene = molecule('C2H4')
        m_ethene = ASE_utils.to_molmod(ethene)
        m_ethene.graph.neighbors
        h_e = huckel.huckel_e([{0, 1}], m_ethene)

        self.assertAlmostEquals(h_e, 2*huckel.alpha + 2*huckel.beta)

