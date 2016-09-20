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
    argparser = ArgumentParser(description="Split file by hit frequency.")
    argparser.add_argument("infile", type = str,  help = "Peptblast output file.")
    argparser.add_argument("head", type=int, help="Size of 'head'")
    argparser.add_argument("tail", type=int, help="Size of 'tail'")
    argparser = argparser.parse_args()

    freqs = defaultdict(int)
    chunks = []

    for ch in grouper(7, open(argparser.infile).readlines()):
        query = ch[1].split('\t')[0]
        freqs[query] += 1
        chunks.append(ch)

    freqs = freqs.items()
    freqs.sort(key = itemgetter(1))

    head = [k for k,v in freqs[:argparser.head]]
    mid  = [k for k,v in freqs[argparser.head:-argparser.tail]]
    tail = [k for k,v in freqs[-argparser.tail:]]

    assert len(head) + len(mid) + len(tail) == len(freqs)

    with open(argparser.infile.replace(".txt", "_head.txt"), "w") as fhead, \
         open(argparser.infile.replace(".txt", "_mid.txt"), "w") as fmid, \
         open(argparser.infile.replace(".txt", "_tail.txt"), "w") as ftail:
        for ch in chunks:
            query = ch[1].split('\t')[0]
            if query in head:
                fhead.writelines(ch)
            elif query in mid:
                fmid.writelines(ch)
            elif query in tail:
                ftail.writelines(ch)
            else:
                assert False, "Something's gone wrong..."

if __name__ == '__main__':
    main()
