#!/usr/bin/python2.7
import numpy as np
from numpy import cos, sin, array, matrix, abs
import os, sys
import pybel
from scipy.spatial import ConvexHull

import utils
import find_existing, list_possible
from prepare import pdb2gmx


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
    # z = ax + by + c
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
    if abs(norm[2] - 1.0) <= 0.001:
        return
    beta = np.arccos(norm[2])
    alpha = np.arcsin(norm[0] / sin(beta))
    # alpha = np.arccos(-norm[1] / sin(beta))
    assert abs(norm[1] + cos(alpha) * sin(beta)) <= 0.001 or abs(norm[1] - cos(alpha) * sin(beta)) <= 0.001
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


def get_line((x0, y0), (x1, y1)):
    return y1 - y0, x0 - x1, (y0 - y1) * x0 + (x1 - x0) * y0


def calc_distance((a, b, c), v):
    return abs(a * v[0] + b * v[1] + c) / np.sqrt(a ** 2 + b ** 2)


def find_the_line(points):
    if len(points) == 2:
        x1, y1 = points[0]
        x2, y2 = points[1]
        k = -(x2 / x1) / (y2 - y1)
        b = (y1 + y2) / 2 - k * (x1 + x2) / 2
        return k, b
    hull = ConvexHull(points)
    best_triple = None
    best_distance = 2007012811.0
    for simplex in hull.simplices:
        longest = 0.0
        longestv = None
        line = get_line(points[simplex][0], points[simplex][1])
        for v in hull.vertices:
            if v in simplex:
                continue
            dist = calc_distance(line, points[v])
            if dist > longest:
                longest, longestv = dist, v
        if longest < best_distance:
            best_distance = longest
            best_triple = (line, longestv)
    assert best_triple is not None
    a, b, c = best_triple[0]
    k = -a / b
    nb = -c / b
    lx, ly = points[best_triple[1]]
    p = a * lx + b * ly + c
    d = best_distance / 2
    if p > 0:
        nb += d * np.sqrt(k ** 2 + 1)
    elif p < 0:
        nb -= d * np.sqrt(k ** 2 + 1)
    else:
        assert False
    return k, nb


def put_residue(res, edge, upper, mirror):
    atom_names = map(lambda x: x[0], list_possible.bond_atoms[res.name][edge])
    atoms = [find_existing.find_atom(res, atom_name) for atom_name in atom_names]
    global testatoms
    testatoms = atoms
    points = array([[_.coords[0], _.coords[1]] for _ in atoms])
    # Should not use linear regression.
    # A = np.vstack([x, np.ones(len(x))]).T
    # k, b = np.linalg.lstsq(A, y)[0]
    k, b = find_the_line(points)

    gamma = np.arctan(k)
    rm = get_rotation_matrix(0, 0, -gamma)
    translate(res, (0, -b, 0))
    rotate(res, rm)
    if (atoms[0].coords[0] > atoms[1].coords[0]) ^ mirror:
        horizontal_flip(res)
    if check_upper(res, atoms) != upper:
        vertical_flip(res)
    mean_y = sum(_.coords[1] for _ in atoms) / len(atoms)
    translate(res, (0, -mean_y + (1.5 if upper else -1.5), 0))


def has_overlay(res1, res2):
    points1, points2 = (array([[_.coords[0], _.coords[1]] for _ in res.atoms]) for res in (res1, res2))
    hull1, hull2 = (ConvexHull(points) for points in (points1, points2))
    for v1 in hull1.vertices:
        p1 = points1[v1]
        l = len(hull2.vertices)
        q = np.cross(points2[hull2.vertices[l - 1]] - p1, points2[hull2.vertices[0]] - p1)
        flag = True
        for i in xrange(l - 1):
            if np.cross(points2[hull2.vertices[i]] - p1, points2[hull2.vertices[i + 1]] - p1) * q < 0:
                flag = False
                break
        if flag:
            return True
    return False


def reflect_point((x0, y0), (k, b)):
    c = (k ** 2 - 1) * x0 + 2 * (b - y0) * k
    x0p = -c / (k ** 2 + 1)
    #    c = abs(k) * (k * x0 + b - y0) * 2 / (k ** 2 + 1)
    #    x0p = x0 + c
    y0p = -(x0p - x0) / k + y0
    return x0p, y0p


def reflect_mutate(res):
    c4p = utils.get_atom_by_name(res, "C4'")
    c5p = utils.get_atom_by_name(res, "C5'")
    to_mutate = map(lambda x: utils.get_atom_by_name(res, x), [
        "H5'1", "H5'2", "O5'", "P", "O1P", "O2P"
    ])
    a, b, c = get_line(c4p.coords[:2], c5p.coords[:2])
    k = -a / b
    bb = -c / b
    for each in to_mutate:
        x, y, z = each.coords
        xp, yp = reflect_point((x, y), (k, bb))
        set_coords(each, (xp, yp, z))


def handle_phosphate(res1, res2):
    if not has_overlay(res1, res2):
        return
    reflect_mutate(res1)
    if not has_overlay(res1, res2):
        return
    reflect_mutate(res2)
    if not has_overlay(res1, res2):
        return
    reflect_mutate(res1)
    if not has_overlay(res1, res2):
        return
    assert False


def try_edge_pair(item, res1, res2):
    bonds = []
    for bond in item["bonds"]:
        atom1 = find_existing.find_atom(res1, bond[0])
        atom2 = find_existing.find_atom(res2, bond[1])
        bonds.append((atom1, atom2))
    move_to_origin(res1)
    move_to_origin(res2)
    rotate_to_xoy(res1)
    rotate_to_xoy(res2)
    iso, edge1, edge2 = item["type"]
    put_residue(res1, edge1, True, False)
    put_residue(res2, edge2, False, True if iso == 't' else False)
    t = -sum(atom1.coords[0] - atom2.coords[0] for atom1, atom2 in bonds) / len(bonds)
    translate(res1, (t, 0, 0))
    p = -sum(atom1.coords[0] + atom2.coords[0] for atom1, atom2 in bonds) / 2 / len(bonds)
    translate(res1, (p, 0, 0))
    translate(res2, (p, 0, 0))
    while len(item["bonds"]) == 2:
        if item["bonds"][0][0] == item["bonds"][1][0]:
            nitro = bonds[0][0]
            res_p = res1
            upper = True
        elif item["bonds"][0][1] == item["bonds"][1][1]:
            nitro = bonds[0][1]
            res_p = res2
            upper = False
        else:
            break
        gamma = None
        x, y, z = nitro.coords
        for neighbour in utils.get_bonded_atoms(nitro):
            if neighbour.GetAtomicNum() == 1:
                continue
            k = (neighbour.GetY() - y) / (neighbour.GetX() - x)
            gamma = np.arctan(k)
            break
        assert gamma is not None
        rm = get_rotation_matrix(0, 0, np.pi / 2 - gamma)
        translate(res_p, (-x, -y, -z))
        rotate(res_p, rm)
        translate(res_p, (x, y, z))
        translate(res_p, (0, (-0.8 if upper else 0.8), 0))
        break


def make_new(item):
    res1c, res2c = item["recs"]
    iso, edge1, edge2 = item["type"]
    ori = str("new_coors/cWW-%c%c.gro" % (res1c, res2c))
    swap = False
    if not os.path.isfile(ori):
        ori = str("new_coors/cWW-%c%c.gro" % (res2c, res1c))
        swap = True
    if not os.path.isfile(ori):
        ori = str("new_coors/tWW-%c%c.gro" % (res1c, res2c))
        swap = False
    assert os.path.isfile(ori)
    moldata = pybel.readfile("gro", ori).next()
    res1, res2 = moldata.residues
    if swap:
        res1, res2 = res2, res1
    assert res1.name == res1c and res2.name == res2c
    new_name = utils.get_name(item)
    filename = "test_new/%s.gro" % new_name

    try_edge_pair(item, res1, res2)
    moldata.localopt()
    moldata.title = new_name
    moldata.write(format="gro", filename=filename, overwrite=True)
    pdb2gmx.editconf(filename, filename, force=True)

    item["nswap"] = swap
    return filename


def main(force=False):
    data_list = list_possible.list_possible()
    for item in data_list:
        if item.get("existing_data") is None and (force or not item.get("created")):
            print utils.get_name(item)
            make_new(item)
            item["created"] = True
    list_possible.save_data(data_list)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        force = sys.argv[1] != "noforce"
    else:
        force = True
    main(force)
