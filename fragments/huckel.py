__author__ = 'clyde'

import numpy as np
from scipy import linalg
from utils import unwind

#params that define the huckel energy
alpha, beta = -11.2, -0.7

def huckel_e(frag, master):
    """Compute the Huckel energy of a fragment"""

    if frag:
        ns = {}
        set_atoms = set(unwind(frag))
        no_atoms = len(set_atoms)

        for ring in frag:
                for i in ring:
                    ns[i] = [n for n in master.graph.neighbors[i] if n in set_atoms]

    else:
        ns = master.graph.neighbors
        no_atoms = len(ns)

    ks = sorted(ns.keys())
    h = np.zeros([len(ns), len(ns)])

    for i, k in enumerate(ks):
        h[i][i] = alpha
        for k2 in ns[k]:
            j = ks.index(k2)
            h[i][j] = beta

    eig_values, eig_vectors = linalg.eigh(h)

    #we doubly occupy all energy levels except if there are an odd number of atoms
    #in which case we singly occupy the homo
    electrons_per_level = []
    last_atom = False
    if no_atoms % 2:
        no_atoms -= 1
        last_atom = True

    for i in range(no_atoms/2):
        electrons_per_level.append(2)
    if last_atom:
        electrons_per_level.append(1)
    for i in range(no_atoms/2):
        electrons_per_level.append(0)

    electrons_per_level = np.array(electrons_per_level)
    return eig_values.dot(electrons_per_level)


def get_atoms_stability(frag, set_atoms, master):
    """Compute the stability of a set of atoms within a fragment"""
    cut_frag = [set(r)-set_atoms for r in frag]
    return huckel_e(cut_frag, master) - huckel_e(frag, master)
