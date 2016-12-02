#!/usr/bin/python2.7
import numpy as np
from numpy import cos, sin, array, matrix
import os
import pybel

import utils
import find_existing


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


def translate(atom, (dx, dy, dz)):
    x, y, z = atom.coords
    set_coords(atom, (x + dx, y + dy, z + dz))


def rotate_to_xoy(res):
    plane = get_plane(res)
    norm = get_norm(plane)
    map(lambda atom: translate(atom, (0, 0, -plane[2])), res.atoms)
    beta = np.arccos(norm[2])
    alpha = np.arccos(-norm[1] / sin(beta))
    assert np.abs(norm[0] - sin(alpha) * sin(beta)) <= 0.001


def try_edge_pair(item, res1, res2):
    bonds = []
    for bond in item["bonds"]:
        atom1 = find_existing.find_atom(res1, bond[0])
        atom2 = find_existing.find_atom(res1, bond[1])
        bonds.append((atom1, atom2))
        for bond in bonds:
            pass
    rotate_to_xoy(res1)
    rotate_to_xoy(res2)


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
    res1, res2 = moldata.residues
    if swap:
        res1, res2 = res2, res1
    assert res1.name == res1c and res2.name == res2c
    new_name = utils.get_name(item)
    filename = "new_coors/%s.gro" % new_name

    try_edge_pair(item, res1, res2)
    moldata.localopt()
    moldata.write(format="gro", filename=filename, overwrite=True)
    return filename
