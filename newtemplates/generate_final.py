#!/usr/bin/python2.7
import os

from newtemplates import find_existing, list_possible
from prepare import pdb2gmx


def find_line(lines, key, cmp=lambda a, b: a == b):
    for i, line in enumerate(lines):
        if cmp(line, key):
            return i
    return -1


def change_num(line, num):
    return line[:6] + "%5d" % num + line[11:]


def handle_existing(item, force):
    if not force and item.get("handled"):
        return
    edata = item["existing_data"]
    if edata["swap"] ^ edata["nswap"] and not item.get("swapped"):
        item["recs"] = item["recs"][1] + item["recs"][0]
        item["type"] = item["type"][0] + item["type"][2] + item["type"][1]
        if item["type"][0] == 't':
            item["bonds"] = list(map(
                lambda x: [x[1], x[0]],
                reversed(item["bonds"])
            ))
            item["mtype"] = item["mtype"][0] + ('' if len(item["mtype"]) == 1 else (
                'a' if item["mtype"][1] == 'b' else 'b'))
        else:
            item["bonds"] = list(map(
                lambda x: [x[1], x[0]],
                item["bonds"]
            ))
            item["mtype"] = str(len(list_possible.bond_atoms[item["recs"][0]][item["type"][1]]) + \
                                len(list_possible.bond_atoms[item["recs"][1]][item["type"][2]]) - 1 - \
                                int(item["mtype"][0]) - 1) + ('' if len(item["mtype"]) == 1 else item["mtype"][1])
        item["swapped"] = True

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
    for i in xrange(atom_start, len(fi_raw) - 2):
        if split_point == -1 and i > atom_start and fi_raw[i - 1][22:26] != fi_raw[i][22:26]:
            split_point = len(fo_raw)
        fo_raw.append(fi_raw[i])
    assert split_point >= 0
    for bond in item["bonds"]:
        fo_raw.append("CONECT %4d %4d\n" % tuple(map(
            lambda raw, atom: int(raw[find_line(raw, atom, lambda a, b: a[12:16].strip() == b)][6:11]),
            (fo_raw, fo_raw[split_point:]), (bond[0], bond[1])
        )))

    with open(new_filename, "w") as fo:
        fo.writelines(fo_raw)
    item["handled"] = True


def main(force=False):
    possible_list = find_existing.main()
    for each in possible_list:
        if each["existing_data"] is not None:
            handle_existing(each, force)
    list_possible.save_data(possible_list)


if __name__ == "__main__":
    main(True)
