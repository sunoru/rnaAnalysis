#!env python2
import shutil
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from non_redundant.utils import *
from non_redundant.repsets import get_nrlist

def collect():
    nrlist_data = get_nrlist('2.140')
    nrlist = nrlist_data['nrlist']
    pwd = os.path.abspath(os.curdir)
    result_dir = os.path.join(DATA_DIR, 'results')
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    os.chdir(TEST_DIR)
    info('Starting collecting the results')
    retested = nrlist_data.get('retested')
    if not retested:
        info('No results to collect')
        return
    diff = open(os.path.join(result_dir, 'diff.txt'), 'w')
    for item in nrlist:
        pdbid = item['rnaid']
        if not retested.get(pdbid):
            continue
        os.chdir(pdbid)
        shutil.copy('result.dat', os.path.join(result_dir, '%s.txt' % pdbid))
        shutil.copy('result_old.dat', os.path.join(result_dir, '%s_old.txt' % pdbid))
        a = b = -1
        with open('diff.txt') as fi:
            for line in fi:
                if line.startswith('<'):
                    a += 1
                elif line.startswith('>'):
                    b += 1
        if a > 0 or b > 0:
            diff.write('%s: +%d -%d\n' % (pdbid, a, b))
        os.chdir('..')
    diff.close()
    info('Collect completed')
    os.chdir(pwd)


if __name__ == '__main__':
    collect()
