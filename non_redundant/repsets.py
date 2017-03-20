#!/usr/bin/python2

from utils import *
import csv

def get_nrlist(version=None, force=False):
    if version is None:
        info("Getting the latest list.")
        x = fetch_raw("http://rna.bgsu.edu/rna3dhub/nrlist/")
        version = findall(r'\">(.*) \(current\)', x)[0]
    datafile = data_file("nrlist-v%s.json" % version)
    if not force and os.path.exists(datafile):
        with open(datafile) as fi:
            return json.load(fi)
    info("Fetching list of representative sets. Version %s" % version)
    csvdata_raw = fetch_raw("http://rna.bgsu.edu/rna3dhub/nrlist/download/%s/2.5A/csv" % version)
    data = {
        "version": version,
        "nrlist": parse_csv(csvdata_raw)
    }
    save_nrlist(data)
    return data


def parse_csv(csvdata_raw):
    alist = []
    for each in csv.reader(csvdata_raw.split('\n')):
        if not each:
            continue
        item = {
            "name": each[0],
            "rnaid": each[1].split('|')[0],
            "representative": each[1],
            "class_members": each[2].split(',')
        }
        alist.append(item)
    return alist


def fetch_rna(rnaid, force=False):
    datadir = os.path.join("testfiles", rnaid)
    if not os.path.isdir(datadir):
        os.mkdir(datadir)
    filename = os.path.join(datadir, "%s.pdb" % rnaid)
    if not force and os.path.exists(filename):
        return filename
    info("Fetching %s pdb file." % rnaid)
    data = fetch_raw("https://files.rcsb.org/download/%s.pdb" % rnaid)
    with open(filename, 'w') as fo:
        fo.write(data)
    return filename


def save_nrlist(nrlist_data):
    datafile = data_file("nrlist-v%s.json" % nrlist_data["version"])
    info("Saving to %s" % datafile)
    with open(datafile, 'w') as fo:
        json.dump(nrlist_data, fo)
    

def fetch_rnas(nrlist, force=False):
    for item in nrlist["nrlist"]:
        if item.get("downloaded"):
            continue
        rnaid = item["rnaid"]
        try:
            fetch_rna(rnaid)
            item["downloaded"] = True
        except e:
            error("Failed to download %s.pdb: %s" % (rnaid, e.msg))
    save_nrlist(nrlist)
