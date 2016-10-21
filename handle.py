#!/usr/bin/python3


import requests
import json

families = ["cWW", "tWW", "cWH", "tWH", "cWS", "tWS", "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"]

def fetch_data():
    for family in families:
        raw = requests.get("http://ndbserver.rutgers.edu/ndbmodule/services/BPCatalog/static/data/%s.json" % family)
        with open("data/%s.json" % family, "wb") as fo:
            fo.write(raw.content)


def load_data():
    data = {"all_types": []}
    for family in families:
        with open("data/%s.json" % family, "r") as fi:
            data[family] = json.load(fi)
            data[family]["iso"] = "cis" if family[0] == 'c' else "trans"
            data[family]["bondtype"] = (family[1], family[2])
        for item in data[family]["items"]:
            data["all_types"].append({
                "flag": False,
                "name": "%s-%s" % (family, item.upper())
            })
    return data
    

def check_item(data, info):
    p1 = ''.join((
        't' if info[2] == "trans" else 'c',
        info[1][0], info[1][2], '-',
        info[0][0], info[0][1],
    ))
    p2 = ''.join((
        't' if info[2] == "trans" else 'c',
        info[1][2], info[1][0], '-',
        info[0][1], info[0][0],
    ))
    p = False
    for i in range(len(data["all_types"])):
        if data["all_types"][i]["name"] == p1 or \
                data["all_types"][i]["name"] == p2:
            data["all_types"][i]["flag"] = True
            p = True
    return p

def main():
    data = load_data()
    with open("RNA-bps.pdb") as fi:
        for line in fi:
            if line.startswith("REMARK"):
                if not check_item(data, line.split()[1:4]):
                    print(line)
    for each in data["all_types"]:
        if not each["flag"]:
            print(each["name"])

# main()

def fetch_coors():
    data = load_data()
    for family in data:
        if family == "all_types":
            continue
        for item_name in data[family]["items"]:
            item = data[family]["items"][item_name]
            if not item["coordinates_exist"]:
                continue
            if item["family"][1].lower() == item["family"][2].lower():
                if item["family"][1] > item["family"][2]:
                    p = '1'
                elif item["family"][1] < item["family"][2]:
                    p = '2'
                else:
                    p = '0'
            else:
                p = '0'
            url = "http://ndbserver.rutgers.edu/ndbmodule/services/BPCatalog/static/data/%s/%s-%c.pdb" % (
                family, item_name.upper(), p
            )
            print(url)
            raw = requests.get(url)
            with open("coors/%s-%s.pdb" % (family, item_name), "wb") as fo:
                fo.write(raw.content)
