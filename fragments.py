import itertools

unwind = lambda l: (e for sl in l for e in sl)

class GraphFrags(object):
    """Class to split a graphitic fragment into subfragments"""
    def __init__(self, master):
        self.master = master
        self.base_frag = tuple([])
        self.shell_lookup = []
        self.ring_lookup = self.get_ring_lookup()
        self.list_rings = self.get_all_rings()

    def set_base_frag(self, frag):
        """Base fragment from which all fragments are generated from and to which all fragments are connected to"""
        self.base_frag = tuple(frag)
        self.shell_lookup = self._get_shell_lookup()

    def length_n_frags(self, initial, n=6):
        """Returns all chains of length n that start at carbon initial"""
        frags = []
        current_frag = initial
        if len(current_frag) >= n:
            return [current_frag]

        neighbor_indices = self.master.graph.neighbors[current_frag[-1]]
        for neighbor_ind in neighbor_indices:
            if neighbor_ind not in current_frag:
                new_frag = current_frag + (neighbor_ind, )
                frags += self.length_n_frags(new_frag, n)
        return frags

    def get_all_rings(self, initial_ring=None, list_rings=None):
        """Recursively generates all rings connected to initial_ring"""
        if not list_rings:
            list_rings = []

        if not initial_ring:
            initial_ring = self.get_rings(0)[0]

        for ring_neighbor in self._get_ring_neighbors(initial_ring):
            if ring_neighbor not in list_rings:
                list_rings.append(ring_neighbor)
                self.get_all_rings(ring_neighbor, list_rings)
        return list_rings

    def get_ring_lookup(self, initial_ring=None):
        """Returns a lookup dictionary specifying the ring neighbors of each ring in the fragment"""
        rings = self.get_all_rings(initial_ring=initial_ring)

        lookup = {}
        for ring in rings:
            lookup[ring] = self._get_ring_neighbors(ring)

        return lookup

    def _get_ring_neighbors(self, initial_ring):
        """"Returns a list of ring neighbors connected to initial_ring"""
        initial_indexes = []
        for i in initial_ring:
            if not any([n in initial_indexes for n in self.master.graph.neighbors[i]]):
                initial_indexes.append(i)

        init_rings = (self.get_rings(i) for i in initial_indexes)

        rings = []
        for r in (r for lr in init_rings for r in lr):
            if r not in rings and r != initial_ring:
                rings.append(r)
        return rings

    def neighbor_combinations(self, neighbors):
        """"Generator of a all possible combinations of the list of rings neighbors"""
        no_neighbors = len(neighbors)

        for i in range(1, no_neighbors+1):
            for n_combination in itertools.combinations(neighbors, i):
                yield n_combination

    def get_rings(self, initial):
        """Returns the list of rings atom initial is part of """
        frags = self.length_n_frags(tuple([initial]))
        rings = []
        for f in frags:
            if initial in self.master.graph.neighbors[f[-1]]:
                ring = frozenset(f)
                if ring not in rings:
                    rings.append(ring)
        return rings

    def redundant_frag(self, frag):
        """Determines whether fragment frag contains holes - i.e. does not include a ring explicitly that is included implicitly because
        all atoms in that ring are included in other rings that are included in the fragment"""
        frag_atoms = set(unwind(frag))

        for ring in self.list_rings:
            if ring not in frag:
                if ring.issubset(frag_atoms):
                    return True
        return False

    def get_frag_neighbors(self, frag):
        """Returns tuple of all neighboring rings for a given fragment"""
        frag_neighbors = []
        for n in unwind(self.ring_lookup[r] for r in frag):
            if n not in frag and n not in frag_neighbors:
                frag_neighbors.append(n)
        return tuple(frag_neighbors)

    def _get_shell_lookup(self):
        """Returns a list of list of rings, element i of this list corresponds to the list of rings for shell i"""
        shell = 1
        allowed_rings = []
        shell_rings = self._get_allowed_rings(shell)
        while shell_rings:
            shell += 1
            allowed_rings.append(shell_rings)
            shell_rings = self._get_allowed_rings(shell)

        return allowed_rings

    def _get_allowed_rings(self, shell, base_frag=None):
        """Gets the nth shell of neighbours from self.base_frag"""

        if not base_frag and not self.base_frag:
            raise(RuntimeError('Require a base fragment from which to compute which rings belong to the shell'))
        elif base_frag:
            shell_base_frag = base_frag
        else:
            shell_base_frag = self.base_frag

        frag_neighbors = []
        for i in range(shell):
            frag_neighbors = self.get_frag_neighbors(shell_base_frag)
            shell_base_frag = shell_base_frag + frag_neighbors

        return frag_neighbors

    def get_allowed_rings(self, shell):
        """Returns the rings belonging to shell shell"""
        try:
            return self.shell_lookup[shell-1]
        except IndexError:
            return tuple([])

    def _gen_fragments(self, base_frag, allowed_rings=None):
        """Generator of fragments build from base_frag and the list of rings specifed in allowed_rings"""
        frag_neighbors = []
        for n in (unwind(self.ring_lookup[r] for r in base_frag)):
            if n not in frag_neighbors:
                frag_neighbors.append(n)

        if allowed_rings is None:
            allowed_rings = frag_neighbors

        allowed_neighbors = [n for n in frag_neighbors if n in allowed_rings]

        for neighbor_combination in self.neighbor_combinations(allowed_neighbors):
            frag = base_frag + neighbor_combination
            if not self.redundant_frag(frag):
                yield frag

    def gen_fragments(self, origin_frag=None, shell=0):
        """Recursively generates fragments starting from origin_frag that are part of shell shell"""
        if not origin_frag:
            origin_frag = self.base_frag
        else:
            origin_frag = tuple(origin_frag)

        if shell == 0:
            yield origin_frag
            shell += 1

        allowed_rings = self.get_allowed_rings(shell)
        generation_1, generation_2 = itertools.tee(self._gen_fragments(origin_frag, allowed_rings))
        for frag in generation_1:
            yield frag

        shell += 1
        for frag in generation_2:
            for next_frag in self.gen_fragments(frag, shell=shell):
                yield next_frag


def get_uniques(frags):
    """Get unique fragments"""
    atoms = lambda e: set(unwind(e))
    list_atoms = []
    for f in frags:
        list_atoms.append(atoms(f))

    unique_atoms = []
    for a in list_atoms:
        if a not in unique_atoms:
            unique_atoms.append(a)

    return unique_atoms


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
    assert len(get_uniques(list_frags)) == len(list_frags)


_test_1a()
_test_1b()