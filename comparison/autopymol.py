# This is a Python 3 script.

from pymol import cmd, stored

RESIDUES = {'A', 'U', 'C', 'G'}

stored.results = []
stored.filename = ""
stored.index = 0

def get_atomnum(tmp):
    assert tmp[0][0] in RESIDUES
    assert tmp[2][0] in RESIDUES
    return tmp[0][1:], tmp[2][1:]

def get_selection(line):
    sel = "resi "
    atomnum = []
    tmp = line.split()
    i = 0
    while tmp[i] != '|':
        i += 1
    if i > 0:
        atomnum.extend(get_atomnum(tmp[:i]))
    atomnum.extend(get_atomnum(tmp[i+1:]))
    sel += '+'.join(set(atomnum))
    return sel

def autocompare(filename=None):
    if filename is not None:
        stored.filename = filename
        stored.index = 0
        with open(filename) as fi:
            stored.results = fi.read().split('\n')[4:]
    while stored.index < len(stored.results):
        line = stored.results[stored.index]
        stored.index += 1
        if line.find("√") >= 0:
            continue
        print(line)
        sel = get_selection(line)
        cmd.select(sel)
        cmd.hide("all")
        cmd.show("lines", sel)
        cmd.label(sel + " and name C1'", '"%s-%s" % (resn,resi)')
        cmd.orient(sel)
        return
    print("Nothing to compare")

cmd.extend("autocompare", autocompare)
cmd.extend("ac", autocompare)
