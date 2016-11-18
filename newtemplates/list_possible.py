#!/usr/bin/python2.7
import json
import os
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


def make_data(iso, edge1, edge2, i, mtype, nuc1, nuc2, bonds):
    return {
        "type": ''.join((iso, edge1, edge2)),
        "mtype": "%s%s" % (i, mtype),
        "recs": nuc1 + nuc2,
        "bonds": bonds
    }


def check(iso, nuc1, nuc2, edge1, edge2):
    result = []
    edges1 = bond_atoms[nuc1][edge1]
    edges2 = bond_atoms[nuc2][edge2] if iso == 'c' else tuple(reversed(bond_atoms[nuc2][edge2]))
    if iso == 't' and nuc1 == 'U' and nuc2 == 'G' and edge1 == 'S' and edge2 == 'W':
        print 0
    i = -len(edges2) + 1
    while i < len(edges1):
        if i >= 0:
            length = min(len(edges1) - i, len(edges2))
        else:
            length = min(len(edges1), i + len(edges2))
        range1 = arange(i, length) if i >= 0 else arange(0, length)
        range2 = arange(0, length) if i >= 0 else arange(-i, length)
        m1, m2 = check_match(edges1, edges2, range1, range2)
        if m1 is None:
            assert m2 is None
            bonds = tuple((edges1[b1][0], edges2[b2][0]) for b1, b2 in izip(range1, range2))
            result.append(make_data(iso, edge1, edge2, i + len(edges2) - 1, '', nuc1, nuc2, bonds))
        else:
            assert m2 is not None
            tmp = []
            for b1, b2 in izip(range1, range2):
                if b1 == m1[0]:
                    assert b2 == m1[1]
                    break
                tmp.append((edges1[b1][0], edges2[b2][0]))
            if tmp:
                bonds = tuple(tmp)
                result.append(make_data(iso, edge1, edge2, i + len(edges2) - 1, 'a', nuc1, nuc2, bonds))
            tmp = []
            for b1, b2 in izip(reversed(range1), reversed(range2)):
                if b1 == m2[0]:
                    assert b2 == m2[1]
                    break
                tmp.append((edges1[b1][0], edges2[b2][0]))
            if tmp:
                bonds = tuple(tmp)
                result.append(make_data(iso, edge1, edge2, i + len(edges2) - 1, 'b', nuc1, nuc2, bonds))
        i += 1

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


datafile = os.path.join(os.path.dirname(__file__), "../data/new_templates.json")


def save_data(result_list):
    with open(datafile, "w") as fo:
        json.dump(result_list, fo)


def list_possible(force=False):
    if not force and os.path.exists(datafile):
        with open(datafile) as fi:
            result_list = json.load(fi)
        return result_list
    result_list = []
    l = len(nuc_names)
    for i in xrange(l):
        for j in xrange(i, l):
            nuc1 = nuc_names[i]
            nuc2 = nuc_names[j]
            result_list.extend(list_pair(nuc1, nuc2))
    save_data(result_list)
    return result_list


if __name__ == "__main__":
    for each in list_possible(True):
        print each
