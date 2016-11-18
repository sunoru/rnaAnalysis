#!/usr/bin/python2.7
import os

from newtemplates import find_existing
from prepare import pdb2gmx


def handle_existing(item):
    edata = item["existing_data"]
    units = [unit.split('|') for unit in edata["units"]]
    pdbfile = "pdbs/%s-%s.pdb" % (units[0][0], edata["id"])
    foname = pdb2gmx.pdb2gro(pdbfile)#, False)
    if not os.path.exists(foname):
        print foname
    else:
        item["handled"] = True


def main():
    possible_list = find_existing.main()
    for each in possible_list:
        if each.get("handled"):
            continue
        if each["existing_data"] is not None:
            handle_existing(each)


if __name__ == "__main__":
    main()
