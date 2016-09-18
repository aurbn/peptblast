#!/usr/bin/python2

import itertools
from operator import itemgetter
from argparse import ArgumentParser
from collections import defaultdict
import re
from sys import exit


def grouper(n, iterable):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk

def desc2wordset(desc, minlen = 3):
    desc = desc.split("OS=")[0] # for hit
    desc = desc.split("[")[0] # for query
    desc = desc.lower()
    words = re.findall(r"[\w]+", desc)
    words = {w for w in words if len(w) >= minlen}
    return words

def main():
    argparser = ArgumentParser(description="Filter by desctiotion")
    argparser.add_argument("infile", type = str,  help = "Peptblast output file.")
    argparser.add_argument('--wordlen', type=int, required=False, default=3, help="Minimal length of a word")
    argparser.add_argument('--minfcommon', type=float, required=False, default=0.5,
                           help="Proportion of common words")
    argparser.add_argument("--out", type=str, required=False, default="filtered.txt",
                           help="Output file.")

    argparser = argparser.parse_args()


    chunks = []
    nch = 0
    for ch in grouper(7, open(argparser.infile).readlines()):
        nch += 1

        try:
            dquery = ch[1].split('\t')[2]
            print dquery
        except IndexError:
            print("No description for %s, leaving" % ch[1].split('\t')[0])
            chunks.append(ch)
            continue

        try:
            dhit = ch[5].split('\t')[2]
        except IndexError:
            print("No description for %s, leaving" % ch[5].split('\t')[0])
            chunks.append(ch)
            continue

        wquery = desc2wordset(dquery)
        whit = desc2wordset(dhit)
        frac = len(wquery.intersection(whit))*1./len(wquery.union(whit))
        if frac <= argparser.minfcommon:
            chunks.append(ch)

    print("%i from %i entries are filtered out" % (nch-len(chunks), nch))

    with open(argparser.out, "w") as out:
        for ch in chunks:
            out.writelines(ch)

if __name__ == "__main__":
    main()