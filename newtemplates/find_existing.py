#!/usr/bin/python2.7

import os

import numpy
import pybel

import list_possible
from prepare import handle, prepare, pdb2gmx


def load_data():
    possible_list = list_possible.list_possible()
    data = handle.load_data()
    data.pop("all_types")
    return possible_list, data


def get_coordfile(item, force=False):
    pdbfile = prepare.fetch_pdb(str(item["pdb"]))
    pdbfile = pdb2gmx.pdb2pdb(pdbfile)
    pdbfile_new = pdbfile.replace("-t", "-%s" % (str(item["id"])))
    if not force and os.path.exists(pdbfile_new):
        return pdbfile_new
    with open(pdbfile) as fi:
        ori = fi.readlines()
    obj_res = []
    for unit in item["units"]:
        unit_data = unit.split('|')
        obj_res.append((unit_data[2], unit_data[4]))
    i = 0
    cchain = cresnum = None
    if obj_res[0][1] == obj_res[1][1]:
        need_modify = "%4d" % (int(obj_res[1][1]) + 1)
    else:
        need_modify = None

    flag = False
    nswap = False

    with open(pdbfile_new, "w") as fo:
        for line in ori:
            chain = line[20:22].strip()
            resnum = line[22:26].strip()
            for ind, each in enumerate(obj_res):
                if chain == each[0] and resnum == each[1]:
                    if not flag:
                        if ind == 1:
                            nswap = True
                        flag = True
                    if resnum != cresnum:
                        cresnum = resnum
                        i = 3
                    if chain != cchain:
                        cchain = chain
                        i = 0
                    if i < 3 and line.find('P') >= 0:
                        i += 1
                        continue
                    if ind == 1 and need_modify:
                        line = line[:22] + need_modify + line[26:]
                    fo.write(line)
                    break
    if need_modify:
        unit_data = item["units"][1].split('|')
        unit_data[4] = need_modify.strip()
        item["units"][1] = '|'.join(unit_data)
    return pdbfile_new, nswap


def find_atom(residue, atomname):
    for atom in residue.atoms:
        if residue.OBResidue.GetAtomID(atom.OBAtom).strip() == atomname:
            return atom
    assert False


def find_res(residues, unit):
    for each in residues:
        if each.OBResidue.GetChain() == unit[0] and each.OBResidue.GetNum() == unit[1]:
            return each
    assert False


def get_distance(coords1, coords2):
    return sum((coords1[i] - coords2[i]) ** 2 for i in xrange(3)) ** 0.5


def get_hydrogen(atom):
    hs = []
    for neighbour in pybel.ob.OBAtomAtomIter(atom):
        bond = atom.GetBond(neighbour)
        if neighbour.GetAtomicNum() == 1 and bond.GetBondOrder() == 1:
            hs.append(neighbour)
    return hs


def get_angle(vec1, vec2, vec3):
    x1, x2, x3 = map(numpy.array, (vec1, vec2, vec3))
    dx1 = x1 - x2
    dx2 = x3 - x2
    w = numpy.cross(dx1, dx2)
    wlen = numpy.linalg.norm(w)
    s = dx1.dot(dx2)
    return numpy.arctan2(wlen, s)


def get_hbond_angle(res1, res2, atom_a, atom_b, bond):
    hs = get_hydrogen(atom_a.OBAtom)
    hs.extend(get_hydrogen(atom_b.OBAtom))
    angle = -1
    for h in hs:
        t = get_angle(atom_a.coords, (h.GetX(), h.GetY(), h.GetZ()), atom_b.coords)
        if t > angle:
            angle = t
    return angle


def calc_score(each, ress, a, b):
    score = 0.0
    for bond in each["bonds"]:
        atom_a = find_atom(ress[a], bond[0])
        atom_b = find_atom(ress[b], bond[1])
        d = get_distance(atom_a.coords, atom_b.coords)
        angle = get_hbond_angle(ress[a], ress[b], atom_a, atom_b, bond)
        score += d - angle
    score /= len(each["bonds"])
    return score


def find_bpdata(possible_list, item, coordfile, nswap):
    print item["id"], "  ",
    family = item["family"].upper()
    family = family[0].lower() + family[1:]
    reversed_family = family[0] + family[2] + family[1]
    sequence = item["sequence"].upper()
    reversed_sequence = sequence[1] + sequence[0]
    units = []
    for unit in item["units"]:
        unit_data = unit.split('|')
        units.append([unit_data[2], int(unit_data[4])])
    assert len(units) == 2

    pdbdata = pybel.readfile("pdb", coordfile).next()
    pdbdata.addh()  # Important. But not saved.
    best_score = 200701281.0
    best_match = None
    for each in possible_list:
        swap = 0
        if each["type"] == family and each["recs"] == sequence:
            swap = 1
        if each["type"] == reversed_family and each["recs"] == reversed_sequence:
            swap |= 2
        if swap == 0:
            continue
        ress = [find_res(pdbdata.residues, units[i]) for i in xrange(2)]

        if swap & 1 == 1:
            score = calc_score(each, ress, 0, 1)
            if score < best_score:
                best_score = score
                best_match = (each, False)
        if swap & 2 == 2:
            score = calc_score(each, ress, 1, 0)
            if score < best_score:
                best_score = score
                best_match = (each, True)

    assert best_match is not None
    if best_match[0].has_key("existing_data"):
        if set(item["units"]) == set(best_match[0]["existing_data"]["units"]):
            return
        if best_score >= best_match[0]["existing_data"]["score"]:
            print item, best_match
            return
        print item, best_match
    best_match[0]["existing_data"] = {
        "id": item["id"],
        "units": item["units"],
        "swap": best_match[1],
        "nswap": nswap,
        "score": best_score
    }
    print best_match[0]["mtype"], best_score


def main(force=False):
    possible_list, data = load_data()
    if not force and possible_list[0].has_key("existing_data"):
        return possible_list
    for family in data:
        for item in data[family]["items"].values():
            if item["pdb"] == "1J8G":
                continue
            if not item["coordinates_exist"] or item["pdb"] == "Modeled":
                continue
            coordfile, nswap = get_coordfile(item)
            find_bpdata(possible_list, item, coordfile, nswap)
    for each in possible_list:
        if not each.has_key("existing_data"):
            each["existing_data"] = None
    list_possible.save_data(possible_list)
    return possible_list


if __name__ == "__main__":
    print main(True)
