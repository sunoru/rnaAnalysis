#!/usr/bin/python2.7
import pybel

def get_name(item):
    return "%s-%s-%s" % (item["type"], item["recs"], item["mtype"])


def atom_name(residue, atom):
    return residue.OBResidue.GetAtomID(atom.OBAtom).strip()

def get_bonded_atoms(atom):
    obatom = atom.OBAtom
    for neighbour in pybel.ob.OBAtomAtomIter(obatom):
        bond = atom.GetBond(neighbour)
        if bond.GetBondOrder() > 0:
            yield neighbour