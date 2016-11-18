#!/usr/bin/python2
# coding=utf-8
import sys

def nt(name):
    try:
        return (name[0], int(name[1:]))
    except:
        return (name, -1)

def input_result(inputfile):
    result = []
    with open(inputfile) as fi:
        for line in fi:
            if not line.startswith("base pair"):
                continue
            tmp = line.split()
            tmp[5]
            t = "%c%c%c" % (
                'c' if tmp[6] == 'cis' else 't',
                tmp[5][0], tmp[5][2]
            )
            result.append((nt(tmp[2]), nt(tmp[4]), t))
    
    result_dssr = []
    for line in sys.stdin:
        if line == '\n':
            continue
        tmp = line.split()
        result_dssr.append((
            nt(tmp[1].split('.')[1]), nt(tmp[2].split('.')[1]),
            tmp[-2]
        ))
    return result, result_dssr

def init_output():
    k = [0, 0, 0]
    op = []
    def tos(r):
        if r is None:
            return ' ' * 16
        return " %s%-3s- %s%-4s%-4s" % (
            r[0][0], r[0][1],
            r[1][0], r[1][1],
            r[2]
        )
    def impl(r1, r2, match=False):
        if r1 is not None: k[0] += 1
        if r2 is not None: k[1] += 1
        if match: k[2] += 1
        op.append("%s%s%s" % (tos(r1), " √ " if match else " | ", tos(r2)))
    def print_output():
        print "There are %d base pairs in our result." % k[0]
        print "There are %d base pairs in DSSR's result." % k[1]
        if k[0] == k[1] == k[2]:
            print "The two results are the same!"
        else:
            print "Only %d lines match." % k[2]
        print "   Our result    |  DSSR's result  "
        for line in op:
            print line
    return impl, print_output

def main(inputfile):
    result, result_dssr = input_result(inputfile)
    output, print_output = init_output()
    i = 0
    j = 0
    while i < len(result) or j < len(result_dssr):
        if i == len(result):
            output(None, result_dssr[j])
            j += 1
            continue
        if j == len(result_dssr):
            output(result[i], None)
            i += 1
            continue
        r1 = result[i]
        r2 = result_dssr[j]
        if r1[0][1] == r2[0][1]:
            if r1[1][1] == r2[1][1]:
                if r1[2] == r2[2]:
                    output(r1, r2, True)
                else:
                    output(r1, r2, False)
                i += 1
                j += 1
                continue
            elif r1[1][1] < r2[1][1]:
                output(r1, None)
                i += 1
                continue
            else:
                output(None, r2)
                j += 1
                continue
        elif r1[0][1] < r2[0][1]:
            output(r1, None)
            i += 1
            continue
        else:
            output(None, r2)
            j += 1
            continue
    print_output()

if __name__ == "__main__":
    assert len(sys.argv) == 2
    main(sys.argv[1])