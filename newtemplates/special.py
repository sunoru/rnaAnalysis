#!/usr/bin/python2.7

from newtemplates.list_possible import list_possible, save_data, bond_atoms, \
    make_data, nuc_names, edge_names


def check_and_add(data_list, iso, nuc1, nuc2, edge1, edge2):
    result = []
    atoms1 = bond_atoms[nuc1][edge1]
    atoms2 = bond_atoms[nuc2][edge2] if iso == 'c' else tuple(reversed(bond_atoms[nuc2][edge2]))
    for i1, atom1 in enumerate(atoms1):
        if atom1[0][0] == 'N' and atom1[1]:
            for i2 in xrange(len(atoms2) - 1):
                if not (atoms2[i2][1] or atoms2[i2 + 1][1]):
                    bonds = ((atom1[0], atoms2[i2][0]), (atom1[0], atoms2[i2 + 1][0]))
                    new_item = make_data(iso, edge1, edge2, i1 + len(atoms2) - i2 - 1.5, 'a', nuc1, nuc2, bonds)
                    existing = None
                    for item in data_list:
                        o = True
                        for p in new_item:
                            if item.get(p) != new_item[p]:
                                o = False
                                break
                        if o:
                            existing = item
                            break
                    if existing:
                        continue
                    data_list.append(new_item)
    for i2, atom2 in enumerate(atoms2):
        if atom2[0][0] == 'N' and atom2[1]:
            for i1 in xrange(len(atoms1) - 1):
                if not (atoms1[i1][1] or atoms1[i1 + 1][1]):
                    bonds = ((atoms1[i1][0], atom2[0]), (atoms1[i1 + 1][0], atom2[0]))
                    new_item = make_data(iso, edge1, edge2, i1 + len(atoms2) - i2 - 0.5, 'b', nuc1, nuc2, bonds)
                    existing = False
                    for item in data_list:
                        o = True
                        for p in new_item:
                            if item.get(p) != new_item[p]:
                                o = False
                                break
                        if o:
                            existing = True
                            break
                    if existing:
                        continue
                    data_list.append(new_item)


def extend_pairs(data_list, nuc1, nuc2):
    for iso in {'t', 'c'}:
        for i1, edge1 in enumerate(edge_names):
            for i2, edge2 in enumerate(edge_names):
                if nuc1 == nuc2 and i2 < i1 or i1 == i2 and nuc_names.index(nuc2) < nuc_names.index(nuc1):
                    continue
                check_and_add(data_list, iso, nuc1, nuc2, edge1, edge2)


def add_special_items(data_list):
    l = len(nuc_names)
    for i in xrange(l):
        for j in xrange(i, l):
            nuc1 = nuc_names[i]
            nuc2 = nuc_names[j]
            extend_pairs(data_list, nuc1, nuc2)


if __name__ == "__main__":
    data_list = list_possible()
    add_special_items(data_list)
    save_data(data_list)
