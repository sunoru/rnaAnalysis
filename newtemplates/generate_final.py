#!/usr/bin/env python2
import os

import itertools

from newtemplates import list_possible, utils
from prepare import make_bps

def mtype_value(x):
    x = float(x)
    if int(x) == x:
        return int(x)
    return x

def add_data(fo, item, coordfile, existing, num):
    itemname = utils.get_name(item)
    print "Working on %s..." % itemname
    with open(coordfile) as fi:
        lines = fi.readlines()
    rmsd = item.get("rmsd", 0.22)

    st = item["mtype"][-1]
    if item["type"][0] == 't':
        if st in {'a', 'b'}:
            reversed_mtype = item["mtype"][:-1] + 'a' if st == 'b' else 'b'
        else:
            reversed_mtype = item["mtype"]
    else:
        if st in {'a', 'b'}:
            reversed_mtype = str(len(list_possible.bond_atoms[item["recs"][0]][item["type"][1]]) +
                                 len(list_possible.bond_atoms[item["recs"][1]][item["type"][2]]) -
                                 2 - mtype_value(item["mtype"][:-1])) + st
        else:
            reversed_mtype = str(len(list_possible.bond_atoms[item["recs"][0]][item["type"][1]]) +
                                 len(list_possible.bond_atoms[item["recs"][1]][item["type"][2]]) -
                                 2 - mtype_value(item["mtype"]))

    fo.write("MODEL%9d\n" % (num + 1))
    fo.write("REMARK {itemname} {mtype}-{reversed_mtype} {rmsd} {artificial}\n".format(
        itemname="%s-%s" % (item["type"], item["recs"]),
        mtype=item["mtype"],
        reversed_mtype=reversed_mtype,
        rmsd=rmsd,
        artificial=0 if existing else 1
    ))
    cryst = tuple(map(lambda x: 10 * float(x), lines[-1].split()))
    fo.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n" % cryst)

    split_point = [j for j in xrange(3, len(lines) - 1) if lines[j - 1].split()[0] != lines[j].split()[0]][0]
    if existing and item["existing_data"]["swap"] or not existing and item.get("nswap"):
        atomlines = list(itertools.chain(lines[split_point:-1], lines[2:split_point]))
    else:
        atomlines = lines[2:-1]
    fr = None
    for i, line in enumerate(atomlines):
        tmp = line.split()
        if fr is None:
            fr = int(tmp[0][:-1])
        coord = map(lambda x: float(x) * 10, tmp[-3:])
        fo.write("ATOM  %5d  %-6s%c%6d    %8.3f%8.3f%8.3f  1.00  0.00\n" % (
            i + 1, tmp[1], tmp[0][-1], 1 if fr == int(tmp[0][:-1]) else 2,
            coord[0], coord[1], coord[2]
        ))

    for bond in item["bonds"]:
        at1, at2 = None, None
        for i, line in enumerate(atomlines):
            if i < split_point - 2:
                if bond[0] == line.split()[1]:
                    assert at1 is None
                    at1 = i + 1
            else:
                if bond[1] == line.split()[1]:
                    assert at2 is None
                    at2 = i + 1
        assert at1 is not None and at2 is not None
        fo.write("CONECT%5d%5d\n" % (at1, at2))

    fo.write("ENDMDL\n")


def main():
    data_list = list_possible.list_possible()
    final_file = open("RNA-bps.pdb", 'w')
    for i, item in enumerate(data_list):
        itemname = utils.get_name(item)
        existing = item.get("existing_data") is not None
        if existing:
            coordfile = item["existing_data"]["coordfile"]
        else:
            coordfile = "test_new/%s.gro" % itemname
        add_data(final_file, item, coordfile, existing, i)
    final_file.close()


if __name__ == "__main__":
    main()
