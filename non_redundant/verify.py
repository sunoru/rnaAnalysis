import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from comparison.auto_compare import analyze
from non_redundant.utils import *
from non_redundant.repsets import get_nrlist, save_nrlist


def verify_one(pdbid):
    os.chdir(pdbid)
    if not os.path.exists('result_old.dat'):
        os.rename('result.dat', 'result_old.dat')
    info('Reanalyzing %s' %pdbid)
    try:
        assert analyze(pdbid) == 0
        os.system('diff result.dat result_old.dat > diff.txt')
        with open('diff.txt') as fi:
            diff = len(fi.readlines())
    except AssertionError as e:
        raise e
    finally:
        os.chdir('..')
    if diff >= 5:
        warn('Different results for %s' % pdbid)
        return False
    return True

def verify_all(force=False):
    nrlist_data = get_nrlist('2.140')
    nrlist = nrlist_data['nrlist']
    if force or 'retested' not in nrlist_data:
        nrlist_data['retested'] = {}
    retested = nrlist_data['retested']
    pwd = os.path.abspath(os.curdir)
    os.chdir(TEST_DIR)
    info('Starting verifying all the pdbs')
    for item in nrlist:
        pdbid = item['rnaid']
        if nrlist_data['tested'].get(pdbid) != 'OK':
            info('%s not tested' % pdbid)
            continue
        if retested.get(pdbid):
            info('Passing %s' % pdbid)
            continue
        try:
            p = verify_one(pdbid)
            retested[pdbid] = 'OK' if p else 'Different'
            save_nrlist(nrlist_data)
        except AssertionError as e:
            why = raw_input('Why? ')
            error(why)
            retested[pdbid] = why
            save_nrlist(nrlist_data)
            cmd = raw_input('Continue? (Y/n) ')
            if cmd.upper() in ['N', 'NO']:
                break
    info('Verify completed')
    os.chdir(pwd)


if __name__ == '__main__':
    verify_all(force=(len(sys.argv) == 2 and sys.argv[1] == 'force'))
