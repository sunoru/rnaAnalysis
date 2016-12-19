#1/usr/bin/python2.7

from newtemplates.make_new import *

data_list = list_possible.list_possible()
# make_new(item)
t = []
for each in data_list:
    # if utils.get_name(each) == "tWW-UG-4":
    if each["existing_data"] is None:
        # make_new(each)
        print utils.get_name(each)
        print each["bonds"]
        t.append(utils.get_name(each))

print "pymol %s" % " ".join(map(lambda x: "%s.gro" % x, t))
