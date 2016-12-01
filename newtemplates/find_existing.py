#!/usr/bin/python2.7
import glob
import os

import numpy
import pybel

import utils
import list_possible
from prepare import handle, prepare, pdb2gmx


def load_data():
    possible_list = list_possible.list_possible()
    coordfiles = glob.glob("coors/*.pdb")
    return possible_list, coordfiles


def add_phosphate(ori, foname):
    with open(ori) as fi1, open(foname) as fi2:
        phosphate = [tuple(map((lambda x: float(x) / 10), fi1.readline()[30:54].split())) for _ in xrange(3)]
        lines = fi2.readlines()
    new_lines = ["%s\n" % os.path.splitext(os.path.split(ori)[-1])[0]]
    new_lines.append("%5d\n" % (int(lines[1]) + 1))
    new_lines.append("%s    P    1%8.3f%8.3f%8.3f\n" % (
        lines[2][:10], phosphate[0][0], phosphate[0][1], phosphate[0][2]
    ))
    new_lines.append("%s  O1P    2%8.3f%8.3f%8.3f\n" % (
        lines[2][:10], phosphate[1][0], phosphate[1][1], phosphate[1][2]
    ))
    new_lines.append("%s  O2P    3%8.3f%8.3f%8.3f\n" % (
        lines[2][:10], phosphate[2][0], phosphate[2][1], phosphate[2][2]
    ))
    i = 3
    for line in lines[2:-1]:
        if line[12:15] in {"H5T", "H3T"}:
            continue
        i += 1
        new_lines.append("%s%3d%s" % (
            line[:17], i, line[20:]
        ))
    new_lines.append(lines[-1])

    with open(foname, "w") as fo:
        fo.writelines(new_lines)
    pdb2gmx.editconf(foname, foname, force=True)

    return foname


def get_coordfile(ori, force=False):
    new_file = ori.replace("coors", "new_coors")
    tmp = os.path.split(ori)[-1]
    new_file = new_file.replace(tmp, tmp[:4] + tmp[4:6].upper() + tmp[6:])

    if not force and os.path.isfile(new_file):
        return new_file
    with open(ori) as fi, open(new_file, "w") as fo:
        for i, line in enumerate(fi):
            if i < 3:
                continue
            fo.write("%sA%s           %c \n" % (line[:21], line[22:-1], line[13]))
    foname = pdb2gmx.pdb2gro(new_file)
    return add_phosphate(ori, foname)


def find_atom(residue, atomname):
    for atom in residue.atoms:
        if utils.atom_name(residue, atom) == atomname:
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


def find_bpdata(possible_list, coordfile):
    print coordfile, "  ",
    family, sequence = os.path.splitext(os.path.split(coordfile)[-1])[0].split('-')
    reversed_family = family[0] + family[2] + family[1]
    reversed_sequence = sequence[1] + sequence[0]

    pdbdata = pybel.readfile("gro", coordfile).next()
    # pdbdata.addh()  # Important. But not saved.
    best_score = 200701281.0
    best_match = None
    ress = pdbdata.residues
    assert ress[0].name + ress[1].name == sequence
    for each in possible_list:
        swap = 0
        if each["type"] == family and each["recs"] == sequence:
            swap = 1
        if each["type"] == reversed_family and each["recs"] == reversed_sequence:
            swap |= 2
        if swap == 0:
            continue

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
        print best_match
    best_match[0]["existing_data"] = {
        "coordfile": coordfile,
        "swap": best_match[1],
        "score": best_score
    }
    print best_match[0]["mtype"], best_score


def main(force=False):
    possible_list, coordfiles = load_data()
    if not force and possible_list[0].has_key("existing_data"):
        return possible_list
    for each in coordfiles:
        coordfile = get_coordfile(each)
        find_bpdata(possible_list, coordfile)

    for each in possible_list:
        if not each.has_key("existing_data"):
            each["existing_data"] = None
    list_possible.save_data(possible_list)
    return possible_list


if __name__ == "__main__":
    print main(True)
