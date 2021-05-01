from sage.all import *
from itertools import tee, combinations

def STSQ(length, cocycle, scomplex):
    """..."""
    answer = set()
    for pair in combinations(cocycle, 2):
        a, b = set(pair[0]), set(pair[1])
        u = a.union(b)
        if len(u) == length:
            u_tuple = tuple(sorted(u))
            if u_tuple in scomplex:
                a_bar, b_bar = a.difference(b), b.difference(a)
                u_bar = sorted(a_bar.union(b_bar))
                index = {}
                for v in a_bar.union(b_bar):
                    pos = u_tuple.index(v)
                    pos_bar = u_bar.index(v)
                    index[v] = (pos + pos_bar) % 2
                index_a = {index[v] for v in a_bar}
                index_b = {index[w] for w in b_bar}
                if index_a == {0} and index_b == {1} or index_a == {1} and index_b == {0}:
                    answer ^= {u_tuple}
    return answer


# noinspection PyProtectedMember
def SQ(self, i):
    """..."""
    P = self.parent()
    scomplex = P.complex()
    if not isinstance(scomplex, SimplicialComplex):
        raise NotImplementedError('This implementation of Steenrod squares is for '
                                  'simplicial complexes only')
    base_ring = P.base_ring()
    if base_ring.characteristic() != 2:
        raise ValueError('Steenrod squares are only defined in characteristic 2')

    ret = P.zero()
    deg_comp = {}
    for index, coeff in self:
        d = deg_comp.get(index[0], {})
        d[index] = coeff
        deg_comp[index[0]] = d

    for j in deg_comp:
        m = j + i
        if not P._graded_indices.get(m, []) or i > j:
            continue
        elt = P._from_dict(deg_comp[j], remove_zeros=False)
        cocycle = [spx[0] for spx in elt.to_cycle()]
        simplices = set((tuple(spx) for spx in scomplex.faces()[m]))
        answer = STSQ(m+1, cocycle, simplices)
        result = {}
        H = scomplex.homology_with_basis(base_ring)
        for gamma_index in H._graded_indices.get(m, []):
            cycle = H._to_cycle_on_basis((m, gamma_index))
            cycle_set = {tuple(spx[0]) for spx in cycle}
            gamma_coeff = len(cycle_set.intersection(answer)) % 2
            if gamma_coeff:
                result[(m, gamma_index)] = gamma_coeff
        ret += P._from_dict(result, remove_zeros=False)
    return ret