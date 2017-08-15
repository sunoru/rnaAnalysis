#!/usr/bin/env python2
import pybel


def get_name(item):
    return str("%s-%s-%s" % (item["type"], item["recs"], item["mtype"]))


def atom_name(residue, atom):
    return residue.OBResidue.GetAtomID(atom.OBAtom).strip()


def get_bonded_atoms(atom):
    obatom = atom.OBAtom
    for neighbour in pybel.ob.OBAtomAtomIter(obatom):
        bond = obatom.GetBond(neighbour)
        if bond.GetBondOrder() > 0:
            yield neighbour


def get_atom_by_name(residue, name):
    name = name.strip()
    for atom in residue.atoms:
        if residue.OBResidue.GetAtomID(atom.OBAtom).strip() == name:
            return atom
    return None
