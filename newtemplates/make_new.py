#!/usr/bin/python2.7
import numpy as np
from numpy import cos, sin, array, matrix
import os
import pybel

import utils
import find_existing, list_possible


def set_coords(atom, (x, y, z)):
    return atom.OBAtom.SetVector(x, y, z)


def get_plane(res):
    a = np.zeros((3, 3))
    b = np.zeros(3)
    n = 0
    for atom in res.atoms:
        name = utils.atom_name(res, atom)
        if name.find('P') >= 0 or name.find("'") >= 0:
            continue
        x, y, z = atom.coords
        a[0][0] += x ** 2
        a[0][1] += x * y
        a[0][2] += x
        a[1][1] += y ** 2
        a[1][2] += y
        b[0] += x * z
        b[1] += y * z
        b[2] += z
        n += 1
    a[1][0] = a[0][1]
    a[2][0] = a[0][2]
    a[2][1] = a[1][2]
    a[2][2] = n
    result = np.linalg.solve(a, b)
    return result


def get_norm((a, b, c)):
    length = np.sqrt(a ** 2 + b ** 2 + 1)
    return (-a / length, -b / length, 1 / length)


def get_rotation_matrix(alpha, beta, gamma):  # in Euler angles
    # I really miss Julia...
    return matrix([
        [cos(alpha) * cos(gamma) - cos(beta) * sin(alpha) * sin(gamma),
         -cos(beta) * cos(gamma) * sin(alpha) - cos(alpha) * sin(gamma), sin(alpha) * sin(beta)],
        [cos(gamma) * sin(alpha) + cos(alpha) * cos(beta) * sin(gamma),
         cos(alpha) * cos(beta) * cos(gamma) - sin(alpha) * sin(gamma), -cos(alpha) * sin(beta)],
        [sin(beta) * sin(gamma), cos(gamma) * sin(beta), cos(beta)]
    ])


def translate_atom(atom, (dx, dy, dz)):
    x, y, z = atom.coords
    set_coords(atom, (x + dx, y + dy, z + dz))


def translate(res, (dx, dy, dz)):
    for atom in res.atoms:
        translate_atom(atom, (dx, dy, dz))


def rotate_atom(atom, rm):
    coords = np.matrix(atom.coords).T
    x, y, z = (rm * coords).T.tolist()[0]
    set_coords(atom, (x, y, z))


def rotate(res, rm):
    for atom in res.atoms:
        rotate_atom(atom, rm)


def rotate_to_xoy(res):
    plane = get_plane(res)
    norm = get_norm(plane)
    translate(res, (0, 0, -plane[2]))
    if np.abs(norm[2] - 1.0) <= 0.001:
        return
    beta = np.arccos(norm[2])
    alpha = np.arcsin(norm[0] / sin(beta))
    # alpha = np.arccos(-norm[1] / sin(beta))
    assert np.abs(norm[1] + cos(alpha) * sin(beta)) <= 0.001
    # assert np.abs(norm[0] - sin(alpha) * sin(beta)) <= 0.001
    rrm = get_rotation_matrix(alpha, beta, 0)
    rm = np.linalg.inv(rrm)
    rotate(res, rm)


def move_to_origin(res):
    center = np.zeros(3)
    for atom in res.atoms:
        center += np.array(atom.coords)
    center /= len(res.atoms)
    translate(res, (-center[0], -center[1], -center[2]))


def flip_atom(atom, direction):
    x, y, z = atom.coords
    if direction == 0:
        # 0 for horizontal
        set_coords(atom, (-x, y, z))
    elif direction == 1:
        set_coords(atom, (x, -y, z))
    else:
        assert False


def vertical_flip(res):
    for atom in res.atoms:
        flip_atom(atom, 1)


def horizontal_flip(res):
    for atom in res.atoms:
        flip_atom(atom, 0)


def check_upper(res, atoms):
    p = 0
    for atom in atoms:
        y = atom.coords[1]
        for neighbour in utils.get_bonded_atoms(atom):
            if neighbour.GetAtomicNum() == 1:
                continue
            if neighbour.GetY() > y:
                p += 1
            else:
                p -= 1
    return p > 0


testpdb = None
iii = 0


def test_show():
    global iii
    iii += 1
    testpdb.write("gro", "test_new/%d.gro" % iii, True)


def put_residue(res, edge, upper, mirror):
    atom_names = map(lambda x: x[0], list_possible.bond_atoms[res.name][edge])
    atoms = [find_existing.find_atom(res, atom_name) for atom_name in atom_names]
    x = [_.coords[0] for _ in atoms]
    y = [_.coords[1] for _ in atoms]
    A = np.vstack([x, np.ones(len(x))]).T
    k, b = np.linalg.lstsq(A, y)[0]
    gamma = -np.arctan(k)
    rm = get_rotation_matrix(0, 0, gamma)
    test_show()
    translate(res, (0, -b, 0))
    test_show()
    rotate(res, rm)
    test_show()
    if not ((atoms[0].coords[0] > atoms[1].coords[0]) ^ mirror):
        horizontal_flip(res)
    if check_upper(res, atoms) != upper:
        vertical_flip(res)
    test_show()
    mean_y = sum(_.coords[1] for _ in atoms) / len(atoms)
    translate(res, (0, -mean_y + (1.5 if upper else -1.5), 0))


def try_edge_pair(item, res1, res2):
    bonds = []
    for bond in item["bonds"]:
        atom1 = find_existing.find_atom(res1, bond[0])
        atom2 = find_existing.find_atom(res2, bond[1])
        bonds.append((atom1, atom2))
    rotate_to_xoy(res1)
    rotate_to_xoy(res2)
    move_to_origin(res1)
    move_to_origin(res2)
    iso, edge1, edge2 = item["type"]
    put_residue(res1, edge1, True, False)
    put_residue(res2, edge2, False, True if iso == 't' else False)
    t = -sum(atom1.coords[0] - atom2.coords[0] for atom1, atom2 in bonds) / len(bonds)
    translate(res1, (t, 0, 0))
    p = -sum(atom1.coords[0] + atom2.coords[0] for atom1, atom2 in bonds) / 2 / len(bonds)
    translate(res1, (p, 0, 0))
    translate(res2, (p, 0, 0))


def make_new(item):
    res1c, res2c = item["recs"]
    iso, edge1, edge2 = item["type"]
    ori = "new_coors/cWW-%c%c.gro" % (res1c, res2c)
    swap = False
    if not os.path.isfile(ori):
        ori = "new_coors/cWW-%c%c.gro" % (res2c, res1c)
        swap = True
    assert os.path.isfile(ori)
    moldata = pybel.readfile("gro", ori).next()
    global testpdb
    testpdb = moldata
    res1, res2 = moldata.residues
    if swap:
        res1, res2 = res2, res1
    assert res1.name == res1c and res2.name == res2c
    new_name = utils.get_name(item)
    # filename = "new_coors/%s.gro" % new_name
    filename = "test_new/%s.gro" % new_name

    try_edge_pair(item, res1, res2)
    # moldata.localopt()
    moldata.write(format="gro", filename=filename, overwrite=True)
    return filename
