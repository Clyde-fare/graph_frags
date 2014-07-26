__author__ = 'clyde'

from utils import get_unique
from fragments import GraphFrags


def _test_0():
    """Tests first shell of neighboring rings correctly identified"""
    from ase.io import read
    import ASE_utils
    from math import pi

    test_mol = read('test1.xyz')
    test_mol.rotate('x', pi/2)
    test_mol.rotate('z', pi/2)

    m_test_mol = ASE_utils.to_molmod(test_mol)
    m_test_mol.set_default_graph()

    g = GraphFrags(m_test_mol)
    g.set_base_frag([g.get_rings(10)[-1]])
    first_shell = g.get_allowed_rings(1)

    assert first_shell == (frozenset({0, 1, 8, 9, 10, 11}),
                           frozenset({9, 10, 17, 18, 19, 20}),
                           frozenset({1, 2, 3, 4, 11, 12}),
                           frozenset({4, 5, 12, 13, 14, 15}),
                           frozenset({19, 20, 21, 22, 25, 26}),
                           frozenset({13, 14, 21, 22, 23, 24}))


def _test_1a():
    """Tests the number of fragments found matches a known case"""
    from ase.io import read
    import ASE_utils
    from math import pi

    test_mol = read('test1.xyz')
    test_mol.rotate('x', pi/2)
    test_mol.rotate('z', pi/2)

    m_test_mol = ASE_utils.to_molmod(test_mol)
    m_test_mol.set_default_graph()

    g = GraphFrags(m_test_mol)
    r = g.get_rings(10)[-1]
    g.set_base_frag([r])
    g_it = g.gen_fragments()

    assert len(list(g_it)) == 112


def _test_1b():
    """Tests no redundant fragments are produced"""
    from ase.io import read
    import ASE_utils
    from math import pi

    base_33_test = read('test1.xyz')
    base_33_test.rotate('x', pi/2)
    base_33_test.rotate('z', pi/2)

    m_33_test = ASE_utils.to_molmod(base_33_test)
    m_33_test.set_default_graph()

    g = GraphFrags(m_33_test)
    r = g.get_rings(10)[-1]
    g.set_base_frag([r])
    g_it = g.gen_fragments()

    list_frags = list(g_it)
    assert len(get_unique(list_frags)) == len(list_frags)

_test_0()
_test_1a()
_test_1b()