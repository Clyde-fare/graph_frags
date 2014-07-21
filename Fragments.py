__author__ = 'clyde'
import itertools
import copy

unwind = lambda l: (e for sl in l for e in sl)

def length_5_frags(master, initial):
    frags = []
    current_frag = initial
    if len(current_frag) >= 6:
        return [current_frag]

    neighbor_inds = master.graph.neighbors[current_frag[-1]]
    for neighbor_ind in neighbor_inds:
        if neighbor_ind not in current_frag:
            new_frag = copy.deepcopy(current_frag)
            new_frag.append(neighbor_ind)
            frags += length_5_frags(master, new_frag)
    return frags

#def in_ring(master, initial_ind):
#    neighbors = [master.graph.neighbors[f[-1]] for f in length_5_frags(master, [initial_ind])]
#    return any([initial_ind in n for n in neighbors])

def get_rings(master, initial):
    frags = length_5_frags(master, [initial])
    rings = []
    for f in frags:
        if initial in master.graph.neighbors[f[-1]] and sorted(f) not in [sorted(r) for r in rings]:
            rings.append(f)
    return rings

def get_ring_neighbors(master, initial_ring):
    #ring_defining_inds = zip(initial_ring, initial_ring[1:])
    init_rings = (get_rings(master, i) for i in initial_ring)

    rings = []
    for r in (r for lr in init_rings for r in lr):
        if sorted(r) not in [sorted(rng) for rng in rings] and sorted(r) != sorted(initial_ring):
            rings.append(r)
    return rings

#def get_all_rings(master, initial_ring=None, list_rings=None):
#    if not list_rings:
#        list_rings = []

#    if not initial_ring:
#        initial_ring = get_rings(master, 0)[0]

#    for ring_neighbor in get_ring_neighbours(master, initial_ring):
#        if sorted(ring_neighbor) not in [sorted(ring) for ring in list_rings]:
#            list_rings.append(ring_neighbor)
#            get_all_rings(master, ring_neighbor, list_rings)
#    return list_rings

#def keep_frag(master, frag_inds):
#    frag = master.graph.get_subgraph(frag_inds)
#    return all([c in in_ring(frag, c) for c in range(len(frag))])

#def extract_fragments(master, current_fragment, list_fragments):

#    for ind in current_fragment:
#        neighbor_inds = master.graph.neighbors[ind]

#        for neighbor_ind in neighbor_inds:
#            if neighbor_ind not in current_fragment:
#                new_frag = current_fragment.copy()
#                new_frag.add(neighbor_ind)

#                if new_frag not in list_fragments:
#                    if keep_frag(master,new_frag):
#                        list_fragments.append(new_frag)

#                    extract_fragments(master, new_frag, list_fragments)

def neighbor_combinations(neighbors):

    list_n_combinations = []
    no_neighbors = len(neighbors)

    for i in range(1, no_neighbors+1):
        list_n_combinations.append(itertools.combinations(neighbors, i))

    return itertools.chain(*list_n_combinations)

def get_ring_fragments(master, list_fragments, initial_frags=None, allowed_rings=None):

    if not initial_frags:
        initial_frags = [[get_rings(master, 0)[0]]]

    if allowed_rings is None:
        neighbors = unwind([get_ring_neighbors(master, r) for initial_frag in initial_frags for r in initial_frag])
        current_rings = [sorted(r) for initial_frag in initial_frags for r in initial_frag]
        allowed_rings = [sorted(r) for r in neighbors if sorted(r) not in current_rings]

    new_frags = []
    for initial_frag in initial_frags:
        new_neighbors = []
        for ring in initial_frag:
            for n in [sorted(rn) for rn in get_ring_neighbors(master, ring) if sorted(rn) in allowed_rings]:
                if n not in new_neighbors:
                    new_neighbors.append(n)

        new_frags += [tuple(initial_frag) + neighbor_combination for neighbor_combination in neighbor_combinations(new_neighbors)]

    if not new_frags:
        return

    allowed_rings = []
    for ring in new_frags[-1]:
        for n in get_ring_neighbors(master, ring):
            sorted_n = sorted(n)
            if sorted_n not in new_frags[-1] and sorted_n not in allowed_rings:
                allowed_rings.append(sorted_n)


    for frag in new_frags:
        list_fragments.append(frag)

    get_ring_fragments(master, list_fragments=list_fragments, initial_frags=new_frags, allowed_rings=allowed_rings)

