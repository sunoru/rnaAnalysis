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
    new_name = "%s-%s-%s" % (item["type"], item["mtype"], item["recs"])
    editconf(foname, "pdbs/%s.pdb" % new_name)
    return
    new_filename = "new_pdbs/%s.pdb" % new_name
    with open(new_filename, "w") as fo:
        fo.write("REMARK {sequence} {bondtype} {rmsd}\n".format(
            sequence=item["recs"],
            bondtype="%s-%s" % (item["type"], item["mtype"]),
            rmsd=0.3
        ))
    item["handled"] = True


def main():
    possible_list = find_existing.main()
    for each in possible_list:
        if each["existing_data"] is not None:
            handle_existing(each)


if __name__ == "__main__":
    main()
