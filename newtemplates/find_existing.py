#!/usr/bin/python2.7

import os
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
    cchain = None
    cresnum = None
    with open(pdbfile_new, "w") as fo:
        for line in ori:
            chain = line[20:22].strip()
            resnum = line[22:26].strip()
            for each in obj_res:
                if chain == each[0] and resnum == each[1]:
                    if resnum != cresnum:
                        cresnum = resnum
                        i = 3
                    if chain != cchain:
                        cchain = chain
                        i = 0
                    if i < 3 and line.find('P') >= 0:
                        i += 1
                        continue
                    fo.write(line)
                    break
    return pdbfile_new


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
    for neighbour in pybel.ob.OBAtomAtomIter(atom):
        bond = atom.GetBond(neighbour)
        if neighbour.GetAtomicNum() == 1 and bond.GetBondOrder() == 1:
            return neighbour
    return None


def check_hbond(res1, res2, atom_a, atom_b, bond):
    h = get_hydrogen(atom_a.OBAtom)



def calc_score(each, ress, a, b):
    score = 0.0
    for bond in each["bonds"]:
        atom_a = find_atom(ress[a], bond[0])
        atom_b = find_atom(ress[b], bond[1])
        d = get_distance(atom_a.coords, atom_b.coords)
        check_hbond(ress[a], ress[b], atom_a, atom_b, bond)
        score += d
    score /= len(each["bonds"])
    return score


def find_bpdata(possible_list, item, coordfile):
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
    best_match[0]["existing_data"] = {
        "id": item["id"],
        "units": item["units"],
        "swap": best_match[1]
    }
    print best_match[0]["mtype"], best_score


def main(force=False):
    possible_list, data = load_data()
    if not force and possible_list[0].has_key("existing_data"):
        return possible_list
    for family in data:
        for item in data[family]["items"].values():
            if not item["coordinates_exist"] or item["pdb"] == "Modeled":
                continue
            coordfile = get_coordfile(item)
            find_bpdata(possible_list, item, coordfile)
    for each in possible_list:
        if not each.has_key("existing_data"):
            each["existing_data"] = None
    list_possible.save_data(possible_list)
    return possible_list


if __name__ == "__main__":
    print main()
