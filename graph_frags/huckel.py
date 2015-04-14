__author__ = 'clyde'

import numpy as np
from scipy import linalg
from utils import unwind

# parameters that define the huckel energy
alpha, beta = -11.2, -0.7


class Huckel(object):
    def __init__(self, master, frag=None):
        #contains the molmod molecule
        self.master = master


        #defines indices of the atoms of interest
        if frag:
            self.frag = frag
        else:
            # set fragment to all sp2 carbons (assumes all atoms with three bonds are sp2 carbons)
            self.frag = [[k for k in self.master.graph.neighbors if len(self.master.graph.neighbors[k]) == 3]]

        self._neighbour_matrix = None
        self._bond_order_matrix = None
        self._electrons_per_level = None
        self._occupied_eigs = None
        self._unoccupied_eigs = None
        self._huckel_eigs = None
        
        self.no_atoms = len(self.neighbour_matrix.keys())
        
    @property
    def neighbour_matrix(self):
        """Get the neighbour matrix associated with the sp2 carbons within a particular fragment"""

        if self._neighbour_matrix is not None:
            return self._neighbour_matrix
        
        neighbour_matrix = {}
        # sp2 atoms in fragment (assumes all atoms with three bonds are sp2 carbons)
        set_atoms = [a for a in set(unwind(self.frag)) if len(self.master.graph.neighbors[a]) == 3]
        
        for a in set_atoms:
            neighbour_matrix[a] = [n for n in self.master.graph.neighbors[a] if n in set_atoms]
        
        self._neighbour_matrix = neighbour_matrix
        return self._neighbour_matrix

    @property
    def bond_order_matrix(self):
        """Convert neighbour matrix into the secular matrix we need to diagonalize"""

        if self._bond_order_matrix is not None:
            return self._bond_order_matrix
            
        list_carbons = sorted(self.neighbour_matrix.keys())
        h = np.zeros([len(self.neighbour_matrix), len(self.neighbour_matrix)])

        for i, k in enumerate(list_carbons):
            h[i][i] = 0
            for k2 in self.neighbour_matrix[k]:
                j = list_carbons.index(k2)
                h[i][j] = 1

        self._bond_order_matrix = h
        return self._bond_order_matrix

    @property
    def huckel_eigs(self):
        """Compute the Huckel eigenvalues/eigenvectors"""

        if self._huckel_eigs is not None:
            return self._huckel_eigs
 
        h = self.bond_order_matrix * beta
        np.fill_diagonal(h, alpha)
        self._huckel_eigs = linalg.eigh(h)
        return self._huckel_eigs

    @property
    def occupied_eigs(self):
        """Return occupied oribtals in order of energy taking
        into account that we doubly occupy orbitals"""
        if self._occupied_eigs is not None:
            return self._occupied_eigs
            
        eigs = self.huckel_eigs[1].T

        alpha_eigs, beta_eigs = [], []

        for ne, eig in zip(self.electrons_per_level, eigs):
            if ne == 2:
                alpha_eigs.append(eig)
                beta_eigs.append(eig)

            elif ne == 1:
                alpha_eigs.append(eig)

        self._occupied_eigs = list(unwind(zip(alpha_eigs, beta_eigs)))
        return self._occupied_eigs
        
    @property
    def unoccupied_eigs(self):
        """Return unoccupied orbitals in order of energy taking
        into account that we doubly occupy orbitals"""
        if self._unoccupied_eigs is not None:
            return self._unoccupied_eigs
            
        eigs = self.huckel_eigs[1].T

        alpha_eigs, beta_eigs = [], []

        for ne, eig in zip(self.electrons_per_level, eigs):
            if ne == 1:
                beta_eigs.append(eig)

            elif ne == 0:
                alpha_eigs.append(eig)
                beta_eigs.append(eig)

        self._unoccupied_eigs = list(unwind(zip(beta_eigs, alpha_eigs)))
        return self._unoccupied_eigs

    @property
    def electrons_per_level(self):
        """Determine number electrons occupying each eigenvector"""
        if self._electrons_per_level is not None:
            return self._electrons_per_level
            
        electrons_per_level = []
        last_atom = False
        if self.no_atoms % 2:
            self.no_atoms -= 1
            last_atom = True

        for i in range(self.no_atoms/2):
            electrons_per_level.append(2)
        if last_atom:
            electrons_per_level.append(1)
        for i in range(self.no_atoms/2):
            electrons_per_level.append(0)

        self._electrons_per_level = np.array(electrons_per_level)
        return self._electrons_per_level

    @property
    def huckel_e(self):
        """Compute the Huckel energy of a fragment"""

        eig_values, eig_vectors = self.huckel_eigs

        # we doubly occupy all energy levels except if there are an odd number of atoms
        # in which case we singly occupy the homo
        electrons_per_level = self.electrons_per_level

        return eig_values.dot(electrons_per_level)

    def bond_order(self, bond_pair):
        """Compute the Huckel bond order associated with a particular bond pair within the fragment"""
        eigen_values, eig_vectors = self.huckel_eigs
        
        list_carbons = sorted(self.neighbour_matrix.keys())
        bond_ids = [list_carbons.index(bond_c) for bond_c in bond_pair]

        bond_order = 0
        for e_n, eig_vec in zip(self.electrons_per_level, eig_vectors.T):
            bond_order += e_n*eig_vec[bond_ids[0]]*eig_vec[bond_ids[1]]

        return bond_order


def get_atoms_stability(frag, set_atoms, master):
    """Compute the stability of a set of atoms within a fragment"""
    cut_frag = [set(r)-set_atoms for r in frag]
    frag_h = Huckel(master, frag)
    cut_frag_h = Huckel(master, cut_frag)
    return cut_frag_h.huckel_e - frag_h.huckel_e

# not trivially invertible i.e. vec_to_frag not easy
def frag_to_vec(frag, master, size=None):
    """Create vector representation of a PAH fragment
    vector is formed of the sorted eigenvalues of the neighbour (bond)  matrix of the fragment"""

    h = Huckel(master, frag).bond_order_matrix()

    no_atoms = len(h)
    if size:
        extras = size - no_atoms
        h = np.pad(h, ((0, extras), (0, extras)), mode='constant', constant_values=0)

    eig_values, eig_vectors = linalg.eigh(h)

    eig_values.sort()
    return eig_values


def frag_to_vec2(frag, master):
    fragment = set(list(unwind(frag)))
    vec = np.zeros(len(master))

    for c in fragment:
        vec[c] = 1

    return vec