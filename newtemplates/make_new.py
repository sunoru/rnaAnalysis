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
    with open(pdbfile) as fi:
        fl = fi.readline()
    foname = pdb2gmx.pdb2gro(pdbfile)  # , False)
    if not os.path.exists(foname):
        print foname
        return
    new_name = "%s-%s-%s" % (item["type"], item["recs"], item["mtype"])
    finame = editconf(foname, "pdbs/%s.pdb" % new_name)
    with open(finame) as fi:
        fi_raw = fi.readlines()
    new_filename = "new_pdbs/%s.pdb" % new_name
    ps = []
    with open("pdbs/%s.pdb" % units[0][0]) as fi:
        tmp = fi.readlines()
    i = find_line(tmp, fl, lambda a, b: a[22:] == b[22:])
    assert i > 0
    while i > 0:
        i -= 1
        if tmp[i][16] == 'B':
            continue
        if tmp[i].startswith("ATOM"):
            ps.append(tmp[i])
            if len(ps) == 3:
                break
    assert ps[2][13] == 'P'
    fo_raw = []
    fo_raw.append("TITLE     %s\n" % new_name)
    fo_raw.append("REMARK {sequence} {bondtype} {rmsd}\n".format(
        sequence=item["recs"],
        bondtype="%s-%s" % (item["type"], item["mtype"]),
        rmsd=0.3
    ))
    fo_raw.append(fi_raw[find_line(fi_raw, "CRYST1", lambda a, b: a.startswith(b))])
    for i in xrange(3):
        fo_raw.append(change_num(ps[2 - i], i + 1))
    atom_start = find_line(fi_raw, "ATOM", lambda a, b: a.startswith(b))
    for i in xrange(atom_start, len(fi_raw) - 2):
        fo_raw.append(change_num(fi_raw[i], i - atom_start + 4))

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
