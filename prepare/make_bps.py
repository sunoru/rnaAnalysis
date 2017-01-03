#!/usr/bin/python2

import os

resorder = {
    'C': [
        ("P", "P"),
        ("OP1", "O1P"),
        ("OP2", "O2P"),
        ("O5'", "O5'"),
        ("C5'", "C5'"),
        ("C4'", "C4'"),
        ("O4'", "O4'"),
        ("C1'", "C1'"),
        ("N1", "N1"),
        ("C6", "C6"),
        ("H6", "H6"),
        ("C5", "C5"),
        ("H5", "H5"),
        ("C4", "C4"),
        ("N4", "N4"),
        ("1H4", "H41"),
        ("2H4", "H42"),
        ("N3", "N3"),
        ("C2", "C2"),
        ("O2", "O2"),
        ("C3'", "C3'"),
        ("C2'", "C2'"),
        ("O2'", "O2'"),
        ("O3'", "O3'")
    ],
    'G': [
        ("P", "P"),
        ("OP1", "O1P"),
        ("OP2", "O2P"),
        ("O5'", "O5'"),
        ("C5'", "C5'"),
        ("C4'", "C4'"),
        ("O4'", "O4'"),
        ("C1'", "C1'"),
        ("N9", "N9"),
        ("C8", "C8"),
        ("H8", "H8"),
        ("N7", "N7"),
        ("C5", "C5"),
        ("C6", "C6"),
        ("O6", "O6"),
        ("N1", "N1"),
        ("H1", "H1"),
        ("C2", "C2"),
        ("N2", "N2"),
        ("1H2", "H21"),
        ("2H2", "H22"),
        ("N3", "N3"),
        ("C4", "C4"),
        ("C3'", "C3'"),
        ("C2'", "C2'"),
        ("O2'", "O2'"),
        ("O3'", "O3'")
    ],
    'A': [
        ("P", "P"),
        ("OP1", "O1P"),
        ("OP2", "O2P"),
        ("O5'", "O5'"),
        ("C5'", "C5'"),
        ("C4'", "C4'"),
        ("O4'", "O4'"),
        ("C1'", "C1'"),
        ("N9", "N9"),
        ("C8", "C8"),
        ("H8", "H8"),
        ("N7", "N7"),
        ("C5", "C5"),
        ("C6", "C6"),
        ("N6", "N6"),
        ("1H6", "H61"),
        ("2H6", "H62"),
        ("N1", "N1"),
        ("C2", "C2"),
        ("H2", "H2"),
        ("N3", "N3"),
        ("C4", "C4"),
        ("C3'", "C3'"),
        ("C2'", "C2'"),
        ("O2'", "O2'"),
        ("O3'", "O3'")
    ],
    'U': [
        ("P", "P"),
        ("OP1", "O1P"),
        ("OP2", "O2P"),
        ("O5'", "O5'"),
        ("C5'", "C5'"),
        ("C4'", "C4'"),
        ("O4'", "O4'"),
        ("C1'", "C1'"),
        ("N1", "N1"),
        ("C6", "C6"),
        ("H6", "H6"),
        ("C5", "C5"),
        ("H5", "H5"),
        ("C4", "C4"),
        ("O4", "O4"),
        ("N3", "N3"),
        ("H3", "H3"),
        ("C2", "C2"),
        ("O2", "O2"),
        ("C3'", "C3'"),
        ("C2'", "C2'"),
        ("O2'", "O2'"),
        ("O3'", "O3'")
    ]
}


def read_residues(lines):
    residues = []
    res = None
    for i, line in enumerate(lines):
        raw = line.split()
        resnum = int(raw[0][:-1])
        if res is None or i > 0 and (resnum != res["resnum"] or raw[0][-1] != res["restype"]):
            res = {
                "resnum": resnum,
                "restype": raw[0][-1],
                "atoms": []
            }
            residues.append(res)
        res["atoms"].append({
            "type": raw[1],
            "num": int(raw[2]),
            "coords": map(lambda x: 10 * float(x), raw[3:])
        })
    if len(residues) == 1:
        atoms = residues[0]["atoms"]
        l = len(atoms)
        assert l & 1 == 0
        residues.append(residues[0].copy())
        residues[0]["atoms"] = atoms[:l / 2]
        residues[1]["atoms"] = atoms[l / 2:]
    return residues


def find_atom(atoms, resod):
    for atom in atoms:
        if atom["type"] == resod[0]:
            return atom, resod[1]
    assert (False)


def find_connect(residues, name, connects_dict={}):
    assert len(residues) == 2
    if not connects_dict:
        with open("hbonds.dat") as fo:
            for line in fo:
                tmp = line.split()
                connects_dict[tmp[0]] = [
                    (tmp[2 * i + 1], tmp[2 * i + 2])
                    for i in xrange(len(tmp) / 2)
                    ]
    bonds = connects_dict[name]
    res1, res2 = residues
    for bond in bonds:
        for atom in res1["atoms"]:
            if atom["type"] == bond[0]:
                atom1 = atom["num"]
                break
        for atom in res2["atoms"]:
            if atom["type"] == bond[1]:
                atom2 = atom["num"]
                break
        yield (atom1, atom2)


def handle_file(fo, filename, i):
    print "Working on %s..." % filename
    family, tmp = filename.split('-')
    sequence = tmp.split('.')[0]
    with open("gros/%s" % filename) as fi:
        lines = fi.readlines()
    if lines[0].startswith('0'):
        rmsd = float(lines[0])
    else:
        rmsd = 0.22

    fo.write("MODEL%9d\n" % (i + 1))
    fo.write("REMARK {sequence} {bondtype} {iso} {rmsd}\n".format(
        sequence=sequence.upper(),
        bondtype=family[1].upper() + family[1].lower() + family[2].upper() + family[2].lower(),
        iso="cis" if family[0] == 'c' else "trans",
        rmsd=rmsd
    ))
    cryst = tuple(map(lambda x: 10 * float(x), lines[-1].split()))
    fo.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n" % cryst)
    residues = read_residues(lines[2:-1])
    q = 0
    resnum = 0
    for res in residues:
        assert (len(res["atoms"]) == len(resorder[res["restype"]]))
        resnum += 1
        res["resnum"] = resnum
        for i in xrange(len(res["atoms"])):
            q += 1
            atom, newname = find_atom(res["atoms"], resorder[res["restype"]][i])
            atom["num"] = q
            atom["type"] = newname
            fo.write("ATOM%7d  %-6s%-6s%-5d%8.3f%8.3f%8.3f  1.00  0.00\n" % (
                q, newname, res["restype"], resnum,
                atom["coords"][0], atom["coords"][1], atom["coords"][2]
            ))
    connects = find_connect(residues, filename.split('.')[0])
    for connect in connects:
        fo.write("CONECT %4d %4d\n" % connect)
    fo.write("ENDMDL\n")


def main():
    coorlist = os.listdir("gros")
    coorlist.sort()

    fo = open("RNA-bps.pdb", "w")

    for i, filename in enumerate(coorlist):
        handle_file(fo, filename, i)

    fo.close()


# with open("t.pdb", "w") as fo:
#    handle_file(fo, "cWH-GG.gro", 0)

if __name__ == "__main__":
    main()
