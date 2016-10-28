#!/usr/bin/python2
# coding=utf-8

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
    if os.path.exists(foname):
        return foname
    with open(finame) as fi, open(foname, 'w') as fo:
        i = 0
        for line in fi:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM"):
                i += 1
                if line.find("HO2'") >= 0:
                    line = line.replace("HO2'", "HO'2")
                if i <= 3 and line.find("P") >= 0:
                    continue
                if line.find("HN3") >= 0:
                    continue # really weird.
                tmp = line.split()
                if len(tmp) > 3 and tmp[3] in rnatypes:
                    fo.write(line)
    return foname

def pdb2gro(pdbname):
    foname = pdbname.replace("pdb", "gro")
    if os.path.exists(foname):
        return
    cmd = "gmx pdb2gmx -f %s -ff amber99sb-ildn -water none -o %s >> /dev/null 2>&1" % (
        pdbname, foname
    )
    print cmd
    os.system(cmd)

def pdb2gmx(pdbid):
    pdbname = pdb2pdb(pdbid)
    pdb2gro(pdbname)

if __name__ == "__pdb2gmx__":
    for each in sys.argv[1:]:
        pdb2gmx(each)
