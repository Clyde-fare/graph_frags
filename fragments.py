import itertools

    
unwind = lambda l: (e for sl in l for e in sl)

class Graph_frags(object):
    def __init__(self, master):
        self.master = master
        self.ring_lookup = self.get_ring_lookup()
        self.list_fragments = []

    def clear(self):
        self.list_fragments = []
        
    def length_6_frags(self, initial):
        frags = []
        current_frag = initial
        if len(current_frag) >= 6:
            return [current_frag]
        
        neighbor_inds = self.master.graph.neighbors[current_frag[-1]]
        for neighbor_ind in neighbor_inds:
            if neighbor_ind not in current_frag:
                #new_frag=copy.copy(current_frag)
                #new_frag.append(neighbor_ind)
                new_frag = current_frag + (neighbor_ind, )
                frags += self.length_6_frags(new_frag)
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
            lookup["-".join(map(str,ring))] = self._get_ring_neighbors(ring)

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
        frags=self.length_6_frags(tuple([initial]))
        rings = []
        for f in frags:
            if initial in self.master.graph.neighbors[f[-1]]:
                ring=sorted(f)
                if ring not in rings:
                    rings.append(ring)
        rings.sort()
        return rings
    
    def get_ring_fragments(self, initial_frags=None, allowed_rings=None):

        if not initial_frags:
            initial_frags = [[self.get_rings(0)[0]]]

        for frag in initial_frags:
            self.list_fragments.append(frag)
            
        if allowed_rings is None:
            neighbors = unwind([self.get_ring_neighbors(r) for initial_frag in initial_frags for r in initial_frag])
            current_rings = [sorted(r) for initial_frag in initial_frags for r in initial_frag]
            allowed_rings = [sorted(r) for r in neighbors if sorted(r) not in current_rings]

        new_frags = []
                        
        for frag in initial_frags:
            new_neighbors = []
            for ring in frag:
                for n in self.get_ring_neighbors(ring):
                    if n in allowed_rings and not n in new_neighbors:
                        new_neighbors.append(n)

            new_frags += [tuple(frag) + neighbor_combination for neighbor_combination in self.neighbor_combinations(new_neighbors)]

        if not new_frags:
            return

        allowed_rings = []
        for ring in new_frags[-1]:
            for n in self.get_ring_neighbors(ring):
                if n not in new_frags[-1] and n not in allowed_rings:
                    allowed_rings.append(n)

        self.get_ring_fragments(initial_frags=new_frags, allowed_rings=allowed_rings)
        
        
def _test_1():
    import ASE_utils
    import math
    from ase.io import read

    base_33_test = read('3_3_test2_flat_graphene.xyz')
    base_33_test.rotate('x', math.pi/2)
    base_33_test.rotate('z', math.pi/2)

    m_33_test = ASE_utils.to_molmod(base_33_test)
    m_33_test.set_default_graph()

    g = Graph_frags(m_33_test)
    r = g.get_rings(10)[-1]
    g.get_ring_fragments(initial_frags=[[r]])
    
    assert len(g.list_fragments) == 112
    
_test_1()