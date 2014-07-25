import itertools

unwind = lambda l: (e for sl in l for e in sl)


class Graph_frags(object):
    def __init__(self, master):
        self.master = master
        self.list_rings = self.get_all_rings()
        self.ring_lookup = self.get_ring_lookup()
        self.list_fragments = []
        self.shell_lookup = []
        self.base_frag = tuple([])

    def clear(self):
        self.list_fragments = []

    def set_base_frag(self, frag):
        self.base_frag = tuple(frag)
        self.shell_lookup = self._get_shell_lookup()

    def length_n_frags(self, initial, n=6):
        frags = []
        current_frag = initial
        if len(current_frag) >= n:
            return [current_frag]

        neighbor_inds = self.master.graph.neighbors[current_frag[-1]]
        for neighbor_ind in neighbor_inds:
            if neighbor_ind not in current_frag:
                new_frag = current_frag + (neighbor_ind, )
                frags += self.length_n_frags(new_frag, n)
        return frags

    def get_all_rings(self, initial_ring=None, list_rings=None):
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
        rings = self.get_all_rings(initial_ring=initial_ring)

        lookup = {}
        for ring in rings:
            lookup["-".join(map(str, ring))] = self._get_ring_neighbors(ring)

        return lookup

    def _get_ring_neighbors(self, initial_ring):
        initial_ring.sort()

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

    def get_ring_neighbors(self, initial_ring):
        return self.ring_lookup["-".join(map(str, initial_ring))]

    def neighbor_combinations(self, neighbors):
        list_n_combinations = []
        no_neighbors = len(neighbors)

        for i in range(1, no_neighbors+1):
            list_n_combinations.append(itertools.combinations(neighbors, i))

        return itertools.chain(*list_n_combinations)

    def get_rings(self, initial):
        frags = self.length_n_frags(tuple([initial]))
        rings = []
        for f in frags:
            if initial in self.master.graph.neighbors[f[-1]]:
                ring = sorted(f)
                if ring not in rings:
                    rings.append(ring)
        rings.sort()
        return rings

    def redundant_frag(self, frag):
        frag_atoms = set(unwind(frag))

        for ring in self.list_rings:
            if ring not in frag:
                if set(ring).issubset(frag_atoms):
                    return True
        return False

    def get_frag_neighbors(self, frag):
        frag_neighbors = []
        for n in unwind(self.get_ring_neighbors(r) for r in frag):
            if n not in frag and n not in frag_neighbors:
                frag_neighbors.append(n)
        return tuple(frag_neighbors)

    def _get_shell_lookup(self):
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
        try:
            return self.shell_lookup[shell-1]
        except IndexError:
            return tuple([])

    def _gen_fragments(self, base_frag, allowed_rings=None):
        frag_neighbors = []
        for n in (unwind([self.get_ring_neighbors(r) for r in base_frag])):
            if n not in frag_neighbors:
                frag_neighbors.append(n)

        if allowed_rings is None:
            allowed_rings = frag_neighbors

        allowed_neighbors = [n for n in frag_neighbors if n in allowed_rings]

        for neighbor_combination in self.neighbor_combinations(allowed_neighbors):
            frag = tuple(base_frag) + neighbor_combination
            if not self.redundant_frag(frag):
                yield frag

    def gen_fragments(self, origin_frag=None, shell=0):

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


def _test_1():
    from ase.io import read
    import ASE_utils
    from math import pi

    base_33_test = read('test1.xyz')
    base_33_test.rotate('x', pi/2)
    base_33_test.rotate('z', pi/2)

    m_33_test = ASE_utils.to_molmod(base_33_test)
    m_33_test.set_default_graph()

    g = Graph_frags(m_33_test)
    r = g.get_rings(10)[-1]
    g.set_base_frag([r])
    g_it = g.gen_fragments()

    assert len(list(g_it)) == 112
    
_test_1()