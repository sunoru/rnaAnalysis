#!/usr/bin/env python2
import pybel
from .utils import *
from .repsets import get_nrlist, save_nrlist


def get_database():
    dbdata = []
    for mol in pybel.readfile('pdb', 'RNA-bps.pdb'):
        remarks = mol.data['REMARK'].split()
        assert len(remarks) == 4
        bondtype = remarks[0]
        subtypes = remarks[1].split('-')
        template = {
            'name': '%s-%s' % (bondtype, subtypes[0]),
            'name_r': '%c%c%c-%c%c-%s' % (bondtype[0], bondtype[2], bondtype[1], bondtype[5], bondtype[4], subtypes[1]),
            'score': float(remarks[2]),
            'modeled': remarks[3] == '1',
            'mol': mol
        }
        dbdata.append(template)
    return dbdata


def mean_all(force=False):
    pass


if __name__ == '__main__':
    mean_all(force=(len(sys.argv) == 1 && sys.argv[1] == 'force'))

