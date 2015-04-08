__author__ = 'clyde'

import numpy as np
from scipy import linalg
from utils import unwind

#params that define the huckel energy
alpha, beta = -11.2, -0.7

def get_neighbour_matrix(frag, master):
    """Get the neighbour matrix associated with sp2 carbons in a particular PAH"""

    if frag:
        ns = {}
        #sp2 atoms in fragment (assuming all atoms are carbons)
        set_atoms = [a for a in set(unwind(frag)) if len(master.graph.neighbors[a]) <= 3]
        for ring in frag:
            for i in ring:
                #include sp2 atoms in the fragment
                if i in set_atoms:
                    ns[i] = [n for n in master.graph.neighbors[i] if n in set_atoms]

    else:
        #dict of atom:[atom_neighbours] for sp2 atoms in total molecule (assumes all atoms are carbons)
        ns = {}

        set_atoms = [k for k in master.graph.neighbors if len(master.graph.neighbors[k]) <= 3]

        for k in set_atoms:
            ns[k] = [n for n in master.graph.neighbors[k] if n in set_atoms]

    ks = sorted(ns.keys())
    h = np.zeros([len(ns), len(ns)])

    for i, k in enumerate(ks):
        h[i][i] = 0
        for k2 in ns[k]:
            j = ks.index(k2)
            h[i][j] = 1

    return h

#not trivially invertible i.e. vec_to_frag not easy
def frag_to_vec(frag, master, size=None):
    """Create vector representation of a PAH fragment
    vector is formed of the sorted eigenvalues of the neighbour (bond)  matrix of the fragment"""
    h = get_neighbour_matrix(frag, master)

    no_atoms = len(h)
    if size:
        extras = size - no_atoms
        h = np.pad(h, ((0,extras),(0,extras)), mode='constant', constant_values=0)

    eig_values, eig_vectors = linalg.eigh(h)

    eig_values.sort()
    return eig_values

def frag_to_vec2(frag, master):
    fragment = set(list(unwind(frag)))
    vec = np.zeros(len(master))

    for c in fragment:
        vec[c] = 1

    return vec

def huckel_e(frag, master):
    """Compute the Huckel energy of a fragment"""

    h = get_neighbour_matrix(frag, master) * beta
    no_atoms = len(h)

    for i in range(no_atoms):
        h[i][i] = alpha

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
