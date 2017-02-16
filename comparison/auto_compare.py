#!/usr/bin/python2
# coding=utf-8
# Usage: python2 auto_compare.py [nohydro] [nophos] [nosugar] [PDBID1] [PDBID2] ...
# If no PDB ID is given, read a PDB ID from each line of stdin.
# It will automatically create a directory named as the PDB ID, download the PDB,
# analyze it with both our tool and DSSR, and finally give the comparison by running `compare.py`.

import codecs
import os
import subprocess
import sys

import requests

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import prepare.pdb2gmx


def writefile(filename):
    return codecs.open(filename, 'w', encoding='UTF-8')

def readfile(filename):
    return codecs.open(filename, 'r', encoding='UTF-8')

def fetch_pdbfile(pdbid):
    filename = "%s.pdb" % (pdbid)
    if os.path.isfile(filename):
        return
    url = "https://files.rcsb.org/download/%s.pdb" % pdbid
    print "Fetching %s" % url
    resp = requests.get(url)
    raw = resp.content.decode("UTF-8")
    with writefile(filename) as fo:
        fo.write(raw)

def editconf(pdbid):
    if os.path.isfile("%s.gro" % pdbid):
        return
    cmd = "gmx editconf -d 1.0 -f %s-t.gro -o %s.gro > /dev/null 2>&1" % (
        pdbid, pdbid
    )
    print cmd
    os.system(cmd)

def analyze(pdbid, hydro=True, phos=True, sugar=True):
    cmd = "gmx rnaAnalysis -f %s.gro -s %s.gro " % (pdbid, pdbid) \
        + "-%shydro -%sphos -%ssugar -g result.dat -o result.xvg > /dev/null 2>&1" % (
            "" if hydro else "no", "" if phos else "no", "" if sugar else "no")
    print cmd
    os.system(cmd)

def fetch_dssr_output(pdbid):
    filename = "dssr.txt"
    if os.path.isfile(filename):
        return
    print "Requesting DSSR output"
    url = "http://dssr.x3dna.org"
    resp = requests.post(url, {"pdbid": pdbid, "submit": "Submit"})
    raw = resp.content.decode("UTF-8")
    with writefile(filename) as fo:
        fo.write(raw)

def get_dssr_result():
    result = []
    with readfile("dssr.txt") as fi:
        for line in fi:
            if line.startswith("      nt1"):
                break
        for line in fi:
            if line == '\n':
                break
            result.append(line)
    return result

def compare(pdbid, options):
    if not os.path.exists(pdbid):
        os.mkdir(pdbid)
    cwd = os.getcwd()
    os.chdir(pdbid)

    fetch_pdbfile(pdbid)
    prepare.pdb2gmx.pdb2gmx(pdbid)
    editconf(pdbid)
    analyze(pdbid, **options)
    fetch_dssr_output(pdbid)
    dssr_result = get_dssr_result()

    print "Comparing %s" % pdbid
    p = subprocess.Popen(["python2", os.path.join("..", os.path.dirname(__file__), "compare.py"), "result.dat"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    for each in dssr_result:
        p.stdin.write(each)
    p.stdin.close()
    with writefile("compare.txt") as fo:
        for line in p.stdout:
            if line.find('√') == -1:
                print line,
            fo.write(line.decode("utf-8"))
    os.chdir(cwd)

def main():
    i = 1
    l = len(sys.argv)
    options = {}
    if l > 1:
        while i < l:
            if sys.argv[i] == "nohydro":
                options["hydro"] = False
            elif sys.argv[i] == "nophos":
                options["phos"] = False
            elif sys.argv[i] == "nosugar":
                options["sugar"] = False
            else:
                break
            i += 1
    if i < l:
        for each in sys.argv[i:]:
            compare(each, options)
    else:
        for line in sys.stdin:
            pdbid = line.strip()
            if pdbid:
                compare(pdbid, options)

main()
