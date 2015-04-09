import unittest
import os

from graph_frags import huckel
from graph_frags.fragments import GraphFrags
from ase_extensions import ase_utils
from ase.io import read
from ase.structure import molecule
from math import pi

__author__ = 'clyde'


class TestHuckel_e(unittest.TestCase):
    def setUp(self):
        #need to use __file__ rather than directly reading test1.xyz because nosetests usually
        #gets run from the directory above tests and would break if we did not tell it that
        #test1.xyz is in the same directory as the python file containing the tests
        test_file1 = os.path.abspath(os.path.dirname(__file__)) + '/test1.xyz'
        test_mol = read(test_file1)
        test_mol.rotate('x', pi/2)
        test_mol.rotate('z', pi/2)

        m_test_mol = ase_utils.to_molmod(test_mol)
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
        m_ethene = ase_utils.to_molmod(ethene)
        m_ethene.graph.neighbors
        ethene_h = huckel.Huckel(m_ethene, [{0,1}])
        h_e = ethene_h.huckel_e()

        self.assertAlmostEquals(h_e, 2*huckel.alpha + 2*huckel.beta)

