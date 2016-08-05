#!/usr/bin/python2

import itertools
from operator import itemgetter
from argparse import ArgumentParser
from collections import defaultdict

def grouper(n, iterable):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk


def main():
    argparser = ArgumentParser(description="Filter and merge enries")
    argparser.add_argument("--1", type = str, dest = "one", required = True,  help = "5+ hits file")
    argparser.add_argument("--2", type = str, dest = "two", required = True,  help = "a/b/c hits file")
    argparser.add_argument("--filter", type = int, required = False,
                           help = "Only use less or equal to cont occurrences")
    argparser.add_argument("--best", action = "store_true", required = False,
                           help = "Select only best match with query/hit combination")
    argparser.add_argument("--out", type=str, required=False, default="merged.txt",
                           help = "Output file.")

    argparser = argparser.parse_args()


    chunks = []
    for ch in grouper(7, open(argparser.one).readlines()):
        freq = int(ch[0].split('\t')[1].split()[1])
        if argparser.filter and freq <= argparser.filter:
            hitid = ch[1].split('\t')[0]
            queryid = ch[5].split('\t')[0]
            metric = int(ch[0].split('\t')[0].split()[0])
            chunks.append((hitid, queryid, metric, ch))


    for ch in grouper(7, open(argparser.two).readlines()):
        freq = int(ch[0].split('\t')[1].split()[1])
        if argparser.filter and freq <= argparser.filter:
            hitid = ch[1].split('\t')[0]
            queryid = ch[5].split('\t')[0]
            a,b,c = ch[0].split('\t')[0].split()[0].split('/')
            metric = int(a) + int(c) - int(b)
            chunks.append((hitid, queryid, metric, ch))

    if argparser.best:
        cdict = defaultdict(list)
        tchunks = []
        for hitid, queryid, metric, ch in chunks:
            cdict[(hitid, queryid)].append((hitid, queryid, metric, ch))
        for k, v in cdict.items():
            maxs = -1000
            maxi = 0
            for i, vi  in enumerate(v):
                if vi[2] > maxs:
                    maxi = i
                    maxs = vi[2]
            tchunks.append(vi)
        chunks = tchunks


    chunks.sort(key = itemgetter(2), reverse=True)

    with open(argparser.out, "w") as out:
        for hitid, queryid, metric, ch in chunks:
            out.writelines(ch)



if __name__ == "__main__":
    main()