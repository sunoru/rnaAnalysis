#!/usr/bin/env python2
import glob
import json
import numpy as np
import pybel
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from non_redundant.utils import *
from non_redundant.repsets import get_nrlist, save_nrlist
from non_redundant.la import jacobi


def get_database():
    info('Reading from the template database')
    dbdata = {}
    for mol in pybel.readfile('pdb', data_file('RNA-bps.pdb')):
        remarks = mol.data['REMARK'].split()
        assert len(remarks) == 4
        bondtype = remarks[0]
        subtypes = remarks[1].split('-')
        template = {
            'name': '%s-%s' % (bondtype, subtypes[0]),
            'name_r': '%c%c%c-%c%c-%s' % (bondtype[0], bondtype[2], bondtype[1], bondtype[5], bondtype[4], subtypes[1]),
            'score': float(remarks[2]),
            'modeled': remarks[3] == '1',
            'mol': mol,
            'connect': [],
            'atomnames': []
        }
        dbdata[template['name']] = template
    with open(data_file('RNA-bps.pdb')) as fi:
        for line in fi:
            if line.startswith('REMARK'):
                tmp = line.split()
                name = '%s-%s' % (tmp[1], tmp[2].split('-')[0])
            elif line.startswith('ATOM'):
                dbdata[name]['atomnames'].append(line[13:17])
            elif line.startswith('CONECT'):
                dbdata[name]['connect'].append(line)
    return dbdata


def get_coords_numatom(mol):
    tmp = mol.data['TITLE'].split()
    score, rmsd = tmp[0], tmp[1]
    swap = mirror = False
    for i in xrange(2, len(tmp)):
        if tmp[i] == 'swap':
            swap = True
        elif tmp[i] == 'mirror':
            mirror = True
        else:
            assert False
    residues = reversed(mol.residues) if swap else mol.residues
    coords = []
    numatom = []
    for residue in residues:
        for atom in residue.atoms:
            x = list(atom.coords)
            if mirror:
                x[2] = -x[2]
            coords.append(x)
        numatom.append(len(residue.atoms))
    return np.array(coords), numatom


def get_coords_mass(mol_or_residue):
    return (
        [atom.coords for atom in mol_or_residue.atoms],
        [atom.atomicmass for atom in mol_or_residue.atoms]
    )


def reset_x(coords, mass):
    l, ndim = coords.shape
    assert ndim == 3
    xcm = mass.dot(coords)
    tm = mass.sum()
    xcm /= tm
    coords -= xcm


def calc_fit_R(coords, tempcoords, mass, sp):
    omega = np.zeros((6, 6))
    om = np.zeros((6, 6))
    d = np.zeros(6)
    # calculate the matrix U
    vh = np.zeros((3, 3))
    vk = np.zeros((3, 3))
    u = np.zeros((3, 3))
    for n in xrange(sp):
        mn = mass[n]
        if mn != 0.0:
            for c in xrange(3):
                xpc = tempcoords[n][c]
                for r in xrange(3):
                    xnr = coords[n][r]
                    u[c][r] += mn * xnr * xpc
    for r in xrange(6):
        for c in xrange(r + 1):
            if r >= 3 and c < 3:
                omega[r][c] = u[r-3][c]
                omega[c][r] = u[r-3][c]
            else:
                omega[r][c] = 0
                omega[c][r] = 0

    irot = jacobi(omega, 6, d, om)

    index = 0
    for j in xrange(2):
        max_d = -1000
        for i in xrange(6):
            if d[i] > max_d:
                max_d = d[i]
                index = i
        d[index] = -10000
        for i in xrange(3):
            vh[j][i] = np.sqrt(2) * om[i][index]
            vk[j][i] = np.sqrt(2) * om[i + 3][index]
    vh[2] = np.cross(vh[0], vh[1])
    vk[2] = np.cross(vk[0], vk[1])

    R = np.zeros((3, 3))
    for r in xrange(3):
        for c in xrange(3):
            for s in xrange(3):
                R[r][c] += vk[s][r]*vh[s][c];
    return R


def do_fit(coords, tempcoords, mass, sp):
    R = calc_fit_R(coords, tempcoords, mass, sp)
    n = coords.shape[0]
    for j in xrange(n):
        vec_old = np.zeros(3)
        for m in xrange(3):
            vec_old[m] = coords[j][m]
        for r in xrange(3):
            coords[j][r] = 0
            for c in xrange(3):
                coords[j][r] += R[r][c] * vec_old[c]


def readpdbfile(pdbfile):
    tmpstr = None
    with open(pdbfile) as fi:
        for line in fi:
            if line.startswith('TITLE'):
                if tmpstr:
                    yield pybel.readstring('pdb', ''.join(tmpstr))
                tmpstr = []
                i = 0
                o = False
            elif line.startswith('ATOM'):
                i += 1
                if i <= 4:
                    chain = str(int(line[25]) + 1)[0]
                    tmpstr.append(line)
                    continue
                if not o and (line[13] == 'P' or line[13:16] == "O5'"):
                    o = True
                if o:
                    line = line[:25] + chain + line[26:]
            tmpstr.append(line)


def mean_pdb(pdbfile, template):
    sumup = [[np.zeros(3), 0] for atom in template['mol'].atoms]
    tempmol = template['mol']
    if 'coords' not in template:
        coords1, mass1 = get_coords_mass(tempmol.residues[0])
        coords2, mass2 = get_coords_mass(tempmol.residues[1])
        template['coords'] = map(np.array, (
            coords1 + coords2,
            coords1[3:] + coords2,
            coords1 + coords2[3:],
            coords1[3:] + coords2[3:]
        ))
        template['mass'] = map(np.array, (
            mass1 + mass2,
            mass1[3:] + mass2,
            mass1 + mass2[3:],
            mass1[3:] + mass2[3:]
        ))
    tempcoordss = template['coords']
    masss = template['mass']

    for mol in readpdbfile(pdbfile):
        coords, numatom = get_coords_numatom(mol)
        t = 0
        if numatom[0] < len(tempmol.residues[0].atoms):
            t += 1
        if numatom[1] < len(tempmol.residues[1].atoms):
            t += 2
        tempcoords = tempcoordss[t]
        mass = masss[t]
        reset_x(coords, mass)
        do_fit(coords, tempcoords, mass, numatom[0])
        for i in xrange(len(coords)):
            if i < numatom[0]:
                j = i + (t & 1) * 3
            else:
                j = i + (t & 1) * 3 + (t & 2) * 3 / 2
            sumup[j][0] += coords[i]
            sumup[j][1] += 1
    return sumup


def mean_one(pdbid, dbdata):
    info('Processing %s' % pdbid)
    os.chdir(pdbid)
    try:
        pdbfiles = [x for x in glob.glob('*.pdb') if x.startswith('cis') or x.startswith('trans')]
        results = {}
        for pdbfile in pdbfiles:
            # Example of `pdbfile`: 'cis_UU-WW-1_UU-WW-3.pdb'
            tmp = pdbfile.split('_')[1].split('-')
            assert len(tmp) == 3
            name = '%c%s-%s-%s' % (pdbfile[0], tmp[1], tmp[0], tmp[2])
            assert name in dbdata
            template = dbdata[name]
            results[name] = mean_pdb(pdbfile, template)
    except AssertionError as e:
        raise e
    finally:
        os.chdir('..')
    return results


def writepdb(fo, dbitem, coords, num=[0]):
    num[0] += 1
    mol = dbitem['mol']
    fo.write('MODEL %8d\n' % num[0])
    fo.write('REMARK %s %s 0.22 1\n' % (
        '-'.join(dbitem['name'].split('-')[:-1]),
        '-'.join((dbitem['name'].split('-')[-1], dbitem['name_r'].split('-')[-1]))
    ))
    fo.write('CRYST1   35.000   35.000   25.000  90.00  90.00  90.00 P 1           1\n')
    for i, (x, y, z) in enumerate(coords):
        fo.write('ATOM  %5d  %4s  %c     2    %8.3f%8.3f%8.3f  1.00  0.00\n' % (
            i + 1,
            dbitem['atomnames'][i],
            'A' if mol.atoms[i].residue.idx == 0 else 'B',
            x, y, z
        ))
    for each in dbitem['connect']:
        fo.write(each)
    fo.write('ENDMDL\n')



def mean_all(force=False):
    nrlist_data = get_nrlist('2.140')
    nrlist = nrlist_data['nrlist']
    dbdata = get_database()
    tmpdb = {}
    pwd = os.path.abspath(os.curdir)
    os.chdir(TEST_DIR)
    info('Starting processing pdbs')
    for item in nrlist:
        pdbid = item['rnaid']
        if nrlist_data['tested'].get(pdbid) != 'OK':
            info('%s not tested' % pdbid)
            continue
        results = mean_one(pdbid, dbdata)
        for name in results:
            dbitem = dbdata[name]
            if name not in tmpdb:
                tmpdb[name] = [[np.array(atom.coords), 1] for atom in dbitem['mol'].atoms]
            for i, (coords, num) in enumerate(results[name]):
                tmpdb[name][i][0] += coords
                tmpdb[name][i][1] += num
    info('Calculating means')
    with open(data_file('new_RNA-bps.pdb'), 'w') as fo:
        for name in dbdata:
            if name not in tmpdb:
                writepdb(fo, dbdata[name], get_coords_mass(dbdata[name]['mol'])[0])
            else:
                coords = [atom[0] / atom[1] for atom in tmpdb[name]]
                writepdb(fo, dbdata[name], coords)
    info('Completed writing new database')
    with open(data_file('tmpdb.json'), 'w') as fo:
        json.dump(tmpdb, fo)


if __name__ == '__main__':
    mean_all(force=(len(sys.argv) == 2 and sys.argv[1] == 'force'))

