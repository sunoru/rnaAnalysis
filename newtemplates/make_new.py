#!/usr/bin/python2.7
import os

from newtemplates import find_existing
from prepare import pdb2gmx


def editconf(finame, foname, print_cmd=True):
    if os.path.exists(foname):
        return foname
    cmd = "gmx editconf -d 1.0 -f %s -o %s > /dev/null 2>&1" % (finame, foname)
    if print_cmd:
        print(cmd)
    os.system(cmd)
    return foname


def find_line(lines, key, cmp=lambda a, b: a == b):
    for i, line in enumerate(lines):
        if cmp(line, key):
            return i
    return -1


def change_num(line, num):
    return line[:6] + "%5d" % num + line[11:]


def handle_existing(item):
    if item.get("handled"):
        return
    edata = item["existing_data"]
    units = [unit.split('|') for unit in edata["units"]]
    pdbfile = "pdbs/%s-%s.pdb" % (units[0][0], edata["id"])
    foname = pdb2gmx.pdb2gro(pdbfile)  # , False)
    if not os.path.exists(foname):
        print foname
        return
    new_name = "%s-%s-%s" % (item["type"], item["recs"], item["mtype"])
    finame = editconf(foname, "pdbs/%s.pdb" % new_name)
    with open(finame) as fi:
        fi_raw = fi.readlines()
    new_filename = "new_pdbs/%s.pdb" % new_name
    if new_name == "cHH-GG-0":
        pass

    fo_raw = []
    fo_raw.append("TITLE     %s\n" % new_name)
    fo_raw.append("REMARK {sequence} {bondtype} {rmsd}\n".format(
        sequence=item["recs"],
        bondtype="%s-%s" % (item["type"], item["mtype"]),
        rmsd=0.3
    ))
    fo_raw.append(fi_raw[find_line(fi_raw, "CRYST1", lambda a, b: a.startswith(b))])
    atom_start = find_line(fi_raw, "ATOM", lambda a, b: a.startswith(b))
    split_point = -1
    for i in xrange(atom_start + 1, len(fi_raw) - 2):
        if fi_raw[i - 1][22:26] != fi_raw[i][22:26]:
            split_point = i
            break
    assert split_point >= 0
    if edata["swap"]:
        atomnr1 = len(fi_raw) - 2 - split_point
        for i in xrange(split_point, len(fi_raw) - 2):
            fo_raw.append(change_num(fi_raw[i], i - split_point + 1))
        for i in xrange(atom_start, split_point):
            fo_raw.append(change_num(fi_raw[i], i - atom_start + atomnr1))
    else:
        atomnr1 = split_point - atom_start
        for i in xrange(atom_start, len(fi_raw) - 2):
            fo_raw.append(fi_raw[i])
    for bond in item["bonds"]:
        fo_raw.append("CONECT %4d %4d\n" % tuple(map(
            lambda raw, atom: int(raw[find_line(raw, atom, lambda a, b: a[12:16].strip() == b)][6:11]),
            (fo_raw, fo_raw[atomnr1:]), (bond[0], bond[1])
        )))

    with open(new_filename, "w") as fo:
        fo.writelines(fo_raw)
    item["handled"] = True


def main():
    possible_list = find_existing.main()
    for each in possible_list:
        if each["existing_data"] is not None:
            handle_existing(each)


if __name__ == "__main__":
    main()
