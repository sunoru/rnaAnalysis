#!/usr/bin/python2.7

from itertools import izip

bond_atoms = {  # No O2' for now.
    'A': {
        'W': (("N6", True), ("N1", False), ("C2", True)),
        'S': (("C2", True), ("N3", False), ("C1'", True)),
        'H': (("N6", True), ("N7", False), ("C8", True))
    },
    'G': {
        'W': (("O6", False), ("N1", True), ("N2", True)),
        'S': (("N2", True), ("N3", False), ("C1'", True)),
        'H': (("O6", False), ("N7", False), ("C8", True))
    },
    'C': {
        'W': (("N4", True), ("N3", False), ("O2", False)),
        'S': (("O2", False), ("C1'", True)),
        'H': (("N4", True), ("C5", True), ("C6", True)),
    },
    'U': {
        'W': (("O4", False), ("N3", True), ("O2", False)),
        'S': (("O2", False), ("C1'", True)),
        'H': (("O4", False), ("C5", True), ("C6", True)),
    }
}

nuc_names = bond_atoms.keys()
edge_names = bond_atoms['A'].keys()


def check_match(edges1, edges2, range1, range2):
    assert len(range1) == len(range2)
    m1 = None
    for i1, i2 in izip(range1, range2):
        if edges1[i1][1] == edges2[i2][1]:
            m1 = (i1, i2)
            break
    m2 = None
    for i1, i2 in izip(reversed(range1), reversed(range2)):
        if edges1[i1][1] == edges2[i2][1]:
            m2 = (i1, i2)
            break
    return m1, m2


def arange(start, length):
    return xrange(start, start + length)


def check(iso, nuc1, nuc2, edge1, edge2):
    result = []
    edges1 = bond_atoms[nuc1][edge1]
    edges2 = bond_atoms[nuc2][edge2] if iso == 't' else reversed(bond_atoms[nuc2][edge2])

    i = -len(edges2) + 1
    while i < len(edges1):
        if i >= 0:
            length = min(len(edges1) - i, len(edges2))
        else:
            length = min(len(edges1), i + len(edges2))
        range1 = arange(i, length) if i >= 0 else arange(i + len(edges2) - 1, length)
        range2 = arange(0, length) if i >= 0 else arange(-i, length)
        m1, m2 = check_match(edges1, edges2, range1, range2)

    return result


def list_pair(nuc1, nuc2):
    result = []
    for iso in {'t', 'c'}:
        for i1, edge1 in enumerate(edge_names):
            for i2, edge2 in enumerate(edge_names):
                if nuc1 == nuc2 and i2 < i1:
                    continue
                result.extend(check(iso, nuc1, nuc2, edge1, edge2))
    return result


def print_result(result_list):
    pass


def list_possible():
    result_list = []
    l = len(nuc_names)
    for i in xrange(l):
        for j in xrange(i, l):
            nuc1 = nuc_names[i]
            nuc2 = nuc_names[j]
            result_list.extend(list_pair(nuc1, nuc2))
    return result_list


def main():
    print_result(list_possible())


if __name__ == "__main__":
    main()
