#!/usr/bin/python2.7

import os
import sys

rnatypes = {"A", "U", "C", "G", "RA", "RU", "RC", "RG", "RA5", "RT5", "RU5", "RC5", "RG5", "RA3", "RT3",
            "RU3", "RC3", "RG3", "RAN", "RTN", "RUN", "RCN", "RGN"}


def pdb2pdb(pdbid):
    if pdbid.endswith(".pdb"):
        finame = pdbid
    else:
        finame = "%s.pdb" % pdbid
    foname = finame.replace(".pdb", "-t.pdb")
    if os.path.isfile(foname):
        return foname
    with open(finame) as fi, open(foname, 'w') as fo:
        cchain = None
        cresnum = None
        starting = True
        for line in fi:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM"):
                if line[16] == 'B':
                    continue
                chain = line[20:22].strip()
                resnum = line[22:26].strip()
                if resnum != cresnum:
                    cresnum = resnum
                    starting = False
                if chain != cchain:
                    cchain = chain
                    starting = True
                if starting and line.find('P') >= 0:
                    continue
                if line.find("HO2'") >= 0:
                    line = line.replace("HO2'", "HO'2")
                if line.find("HN3") >= 0:
                    continue  # really weird.
                tmp = line.split()
                if len(tmp) > 3 and tmp[3] in rnatypes:
                    fo.write(line)
    return foname


def pdb2gro(pdbname, print_cmd=True):
    foname = pdbname[:-3] + "gro"
    if os.path.isfile(foname):
        return foname, 0
    cmd = "gmx pdb2gmx -f %s -ff amber99sb-ildn -water none -p %s -i %s -o %s >> /dev/null 2>&1" % (
        pdbname, pdbname[:-3] + "top", pdbname[:-3] + "itp", foname
    )
    if print_cmd:
        print cmd
    ret = os.system(cmd)
    return foname, ret


def pdb2gmx(pdbid):
    pdbname = pdb2pdb(pdbid)
    return pdb2gro(pdbname)[1]


def editconf(finame, foname, force=False, print_cmd=True):
    if not force and os.path.isfile(foname):
        return foname
    cmd = "gmx editconf -d 1.0 -f %s -o %s > /dev/null 2>&1" % (finame, foname)
    if print_cmd:
        print(cmd)
    os.system(cmd)
    return foname


if __name__ == "__main__":
    for each in sys.argv[1:]:
        pdb2gmx(each)
