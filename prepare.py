#!/usr/bin/python3

import glob
import json
import os
import pickle
import requests

def hasitem(data, item):
    if item["res1"] is None:
        return False
    for each in data.values():
        if each["pdb"] == item["pdb"] and \
                {each["res1i"], each["res2i"]} == {item["res1i"], item["res2i"]}:
            return True
    return False

def usable_data(ori, force_fetch=False):
    datapath = "data/data.json"
    redpath = "data/redundant.dat"
    if not force_fetch and os.path.exists(datapath):
        with open(datapath) as fi:
            rawdata = fi.read()
        with open(redpath, "rb") as fi:
            rawred = fi.read()
        return json.loads(rawdata), pickle.loads(rawred)
    import handle
    ori = handle.load_data()
    data = dict()
    redundant = []
    for family in ori:
        if family == "all_types":
            continue
        for itemname in ori[family]["items"]:
            item = ori[family]["items"][itemname]
            if not item["coordinates_exist"]:
                continue
            name = "%s-%s" % (family, itemname.upper())
            if item["pdb"] == "Modeled":
                res1 = res1i = res2 = res2i = None
            else:
                res1, res1i = item["units"][0].split('|')[3:5]
                res1i = int(res1i)
                res2, res2i = item["units"][1].split('|')[3:5]
                res2i = int(res2i)
            newitem = {
                "pdb": item["pdb"],
                "family": family,
                "sequence": itemname.upper(),
                "res1": res1,
                "res1i": res1i,
                "res2": res2,
                "res2i": res2i
            }
            if hasitem(data, newitem):
                redundant.append("%s-%s" % (family, itemname))
                continue
            data[name] = newitem
    with open(datapath, "w") as fo:
        json.dump(data, fo)
    with open(redpath, "wb") as fo:
        pickle.dump(redundant, fo)
    return data, redundant

# No use.
def fetch_pdbs(data):
    for name in data:
        print(name, end="...")
        item = data[name]
        pdbfile = "pdbs/%s.pdb" % item["pdb"]
        if os.path.exists(pdbfile):
            print("exists")
            continue
        resp = requests.get("https://files.rcsb.org/download/%s.pdb" % item["pdb"])
        raw = resp.content.decode("UTF-8")
        with open(pdbfile, "w") as fo:
            fo.write(raw)
        print("done")

def remove_redundant(redundant):
    for each in glob.glob("coors/*.pdb"):
        if each[6:-4] in redundant:
            os.rename("coors/%s" % each, "coors/%s.bak" % each)

def editconf():
    for each in glob.glob("coors/*.pdb"):
        cmd = "gmx editconf -d 1.0 -f %s -o gros/%s > /dev/null 2>&1" % (each, each[6:-4])
        print(cmd)
        os.system(cmd)

# data, redundant = usable_data()
# remove_redundant(redundant)
editconf()
