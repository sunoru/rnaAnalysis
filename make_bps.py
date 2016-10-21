#!/usr/bin/python2

import os

resorder = {
    'C': (
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
    ),
    'G': (
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
    ),
    'A': (
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
    ),
    'U': (
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
        ("O3'", "O3'"),
    )
}

coorlist = os.listdir("gros")

fo = open("RNA-bps.pdb", "w")

for i in range(len(coorlist)):
    filename = coorlist[i]
    family, tmp = filename.split('-')
    sequence = tmp.split('.')[0]
    with open("gros/%s" % coorlist[i]) as fi:
        lines = fi.readlines()

    fo.write("MODEL%9d\n", i + 1)
    fo.write("REMARK {sequence} {bondtype} {iso} {rmsd}".format(
        sequence = sequence.upper(),
        bondtype = family[1].upper() + family[1].lower() + family[2].upper() + family[2].lower(),
        iso = "cis" if family[0] == 'c' else "trans",
        rmsd = 0.3
    ))
    cryst = tuple(map(lambda x: 10 * float(x), lines[-1].split()))
    fo.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1" % cryst)

fo.close()
