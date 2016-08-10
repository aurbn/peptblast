#!/usr/bin/python2

import itertools
from operator import itemgetter
from argparse import ArgumentParser
from collections import defaultdict
from shutil import copy
from os import unlink


def grouper(n, iterable):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk


def main():
    argparser = ArgumentParser(description="Count id accurencies")
    argparser.add_argument("infile", type = str,  help = "Peptblast output file.")
    argparser.add_argument("--count", type = str, choices = ["query", "hit"], default="query",
                           help = "Id to count")
    argparser.add_argument("--nosort", action="store_true", help = "Doesnt sort an d modify input file.")
    argparser.add_argument("--table", type=str, required=False, default="table.txt",
                           help="Output table of occurrences.")
    argparser.add_argument('--nobk', action='store_true', required=False, help='Do not save backup file.')
    argparser.add_argument("--organisms", type=str, required=False,
                           help="Output table for organism frequencies.")
    argparser.add_argument("--le", type=int, required=False,
                           help="Print only hits occurred less or equal times.")
    argparser = argparser.parse_args()


    orgfreqs = defaultdict(int)
    freqs = defaultdict(int)
    descrs = {}
    chunks = defaultdict(list)
    for ch in grouper(7, open(argparser.infile).readlines()):
        if argparser.count == "query":
            id_ = ch[1].split('\t')[0]
        elif argparser.count == "hit":
            id_ = ch[5].split('\t')[0]
        else:
            exit(1)
        if argparser.organisms:
            org = ch[5].split('\t')[2]
            if '[' in org:
                org = org.split('[')[1]
                org.replace(']', "")
                orgfreqs[org] += 1

        freqs[id_] += 1
        descrs[id_] = ch[1].split('\t')[2].replace('\n', "")
        chunks[id_].append(ch)

    freqs = freqs.items()
    freqs.sort(key=itemgetter(1))

    if argparser.table:
        with open(argparser.table, "w") as tab:
            acc = 0
            tab.write("Count\tCountSum\tID\tDescription\n")
            for id_, f in freqs:
                acc += f
                tab.write("%i\t%i\t%s\t%s\n" % (f, acc, id_, descrs[id_]))

    if argparser.organisms:
        orgfreqs = orgfreqs.items()
        orgfreqs.sort(key = itemgetter(1), reverse=True)
        with open(argparser.organisms, "w") as oout:
            oout.write("Orgamism\tCount\n")
            for fq in orgfreqs:
                oout.write("%s\t%i\n" % fq)

    if not argparser.nosort:
        if not argparser.nobk:
            copy(argparser.infile, argparser.infile+".bk")
        unlink(argparser.infile)
        freqs = [f for f in freqs if (not argparser.le) or f[0] <= argparser.le]

        with open(argparser.infile, "w") as out:
            for id_, count in freqs:
                for chunk in chunks[id_]:
                    out.writelines(chunk[0].replace('\n', '\t') + "occured %i times\n" % count)
                    out.writelines(chunk[1:])


if __name__ == '__main__':
    main()