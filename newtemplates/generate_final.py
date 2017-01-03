import os

from newtemplates import list_possible, utils


def add_data(fo, item, coordfile, existing, i):
    itemname = utils.get_name(item)
    print "Working on %s..." % itemname
    with open(coordfile) as fi:
        lines = fi.readlines()
    rmsd = item.get("rmsd", 0.22)

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
    data_list = list_possible.list_possible()
    final_file = open("RNA-bps.pdb", 'w')
    for i, item in enumerate(data_list):
        itemname = utils.get_name(item)
        existing = item["existing_data"] is not None
        if existing:
            coordfile = "new_coors/%s.gro" % '-'.join(itemname.split('-')[-1])
        else:
            coordfile = "test_new/%s.gro" % itemname
        add_data(final_file, item, coordfile, existing, i)
    final_file.close()


if __name__ == "__main__":
    main()
