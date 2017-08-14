#/usr/bin/python2

from .utils import *
from .repsets import get_nrlist, save_nrlist
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from comparison.auto_compare import editconf, analyze
from prepare.pdb2gmx import pdb2gmx


def testone(pdbid):
    info('Testing %s' % pdbid)
    os.chdir(pdbid)
    try:
        assert pdb2gmx(pdbid) == 0
        assert editconf(pdbid) == 0
        assert analyze(pdbid) == 0
    except AssertionError as e:
        raise e
    finally:
        os.chdir('..')


def testall():
    nrlist_data = get_nrlist('2.140')
    nrlist = nrlist_data['nrlist']
    pwd = os.path.abspath(os.curdir)
    os.chdir(TEST_DIR)
    info('Starting testall')
    for item in nrlist:
        pdbid = item['rnaid']
        if item.get('tested'):
            info('Passing %s' % pdbid)
            continue
        if not item.get("downloaded"):
            info('%s not downloaded' % pdbid)
            continue
        try:
            testone(pdbid)
            item['tested'] = True
        except AssertionError as e:
            if raw_input('Continue? (Y/n)').upper() == 'N':
                break
        finally:
            save_nrlist(nrlist_data)
    os.chdir(pwd)
