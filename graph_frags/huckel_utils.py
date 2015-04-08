__author__ = 'clyde'

from itertools import islice, chain
from multiprocessing.pool import Pool
from itertools import izip
from graph_frags.fragments import GraphFrags
from graph_frags.huckel import get_atoms_stability
from ase_extensions.ase_utils import to_molmod

import numpy as np

def batch(iterable, size):
    sourceiter = iter(iterable)
    while True:
        batchiter = islice(sourceiter, size)
        yield chain([batchiter.next()], batchiter)

def parallel_score_it(chunk_it, score_f, ncpus=2):
    p = Pool(ncpus)
    for chunk in chunk_it:
        for score in p.map(score_f, chunk):
            yield score
    p.close()

def frag_to_ase(frag, a_master):
    """Converts a list of ring fragments to an ase object"""
    #contiguous fragments will generate sensible looking PAHs, neighbouring atoms in the master that are not in the fragment become fragment
    #H atoms this means fragments with holes end up looking weird with unphysical H atoms embedded into them.
    m_master = to_molmod(a_master)
    frag_a_is = {a for r in frag for a in r}
    frag_link_is = {n for i in frag_a_is for n in m_master.graph.neighbors[i] if n not in frag_a_is}

    ase_link_atoms = a_master[list(frag_link_is)]
    ase_core_atoms = a_master[list(frag_a_is)]

    ase_link_atoms.set_chemical_symbols(['H' for i in range(len(ase_link_atoms))])
    ase_frag = ase_core_atoms + ase_link_atoms
    return ase_frag

def get_neighbor_pairs(frag):
    return {frozenset({i,n}) for i in frag.graph.neighbors for n in frag.graph.neighbors[i] if frag.numbers[i]==6 and frag.numbers[n]==6}


from itertools import tee
unwind = lambda ll: [e for l in ll for e in l]

#this doesn't work because p.map doesn't work on methods!
class Frag_data(object):
    def __init__(self, mol, central_pair):

        self.mol=mol
        self.central_pair = central_pair
        self.g = GraphFrags(self.mol)
        self.base_r = next([r for r in self.g.get_all_rings() if self.central_pair in r])
        self.full_pair_e = get_atoms_stability(self.g.get_all_rings(), set_atoms=central_pair, master=self.g.master)

        self.g.set_base_frag([self.base_r])
        self.frags1, self.frags2 = tee(self.g.gen_fragments())

        chunk_it = batch(self.frags1, 2000)

        self.master_no_atoms = len(mol.numbers)
        self.pair_es = self.parallel_pair_e_it(chunk_it)

    def mod_pairs(self, frag):
        return get_atoms_stability(frag, set_atoms=self.central_pair, master=self.g.master)

    def parallel_pair_e_it(self, chunk_it):
        p = Pool(4)
        for chunk in chunk_it:
            for pair_e in p.map(self.mod_pairs, chunk):
                yield pair_e
        p.close()

    def get_stats(self):
        best = [[(0,((0,),((0,))))] for k in range(1,self.master_no_atoms+1)]
        worst = [[(1,((0,),((0,))))] for k in range(1,self.master_no_atoms+1)]

        max_length = 100
        n_frags = 0

        for pair_e, frag in izip(self.pair_es, self.frags2):
            n_frags+=1
            ac = 1/(1+abs(pair_e - self.full_pair_e))
            no_c = len(frozenset(unwind(frag)))

            if ac > best[no_c][-1][0]:
                if len(best[no_c]) >= max_length:
                    best[no_c].pop()
                best[no_c].append((ac, frag))
                best[no_c].sort(key=lambda e: e[0],reverse=True)

            if ac < worst[no_c][-1][0]:
                if len(worst[no_c]) >= max_length:
                    worst[no_c].pop()
                worst[no_c].append((ac, frag))
                worst[no_c].sort(key=lambda e: e[0])

        return n_frags, best, worst