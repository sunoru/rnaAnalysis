#!/usr/bin/env python2
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from comparison.auto_compare import editconf, analyze
from prepare.pdb2gmx import pdb2gmx
from non_redundant.utils import *
from non_redundant.repsets import get_nrlist, save_nrlist


def testone(pdbid, force=False):
    info('Testing %s' % pdbid)
    os.chdir(pdbid)
    try:
        assert pdb2gmx(pdbid, force=force) == 0
        assert editconf(pdbid, force=force) == 0
        assert analyze(pdbid, pdbFiles=True) == 0
    except AssertionError as e:
        raise e
    finally:
        os.chdir('..')


def testall(force=False):
    nrlist_data = get_nrlist('2.140')
    nrlist = nrlist_data['nrlist']
    if force or 'tested' not in nrlist_data:
        nrlist_data['tested'] = {}
    tested = nrlist_data['tested']
    downloaded = nrlist_data['downloaded']
    pwd = os.path.abspath(os.curdir)
    os.chdir(TEST_DIR)
    info('Starting testall')
    for item in nrlist:
        pdbid = item['rnaid']
        if tested.get(pdbid):
            info('Passing %s' % pdbid)
            continue
        if not downloaded.get(pdbid):
            info('%s not downloaded' % pdbid)
            continue
        try:
            testone(pdbid, force)
            tested[pdbid] = 'OK'
            save_nrlist(nrlist_data)
        except AssertionError as e:
            why = raw_input('Why? ')
            error(why)
            tested[pdbid] = why
            save_nrlist(nrlist_data)
            cmd = raw_input('Continue? (Y/n) ')
            if cmd.upper() in ['N', 'NO']:
                break
    os.chdir(pwd)


if __name__ == '__main__':
    testall(force=(len(sys.argv) == 2 and sys.argv[1] == 'force'))
