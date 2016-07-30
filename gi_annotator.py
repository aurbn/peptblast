#!/usr/bin/python2
import sys
import os
import tempfile
import shutil
from argparse import ArgumentParser
import itertools
from Bio import Entrez


def grouper(n, iterable):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk


def main():
    argparser = ArgumentParser(description="Add gi annotation to peptblast output")
    argparser.add_argument("infile", type = str,  help = "Peptblast output file.")
    argparser.add_argument('--nobk', action='store_true', required=False, help='Do not save backup file.')
    argparser = argparser.parse_args()

    outfile = tempfile.NamedTemporaryFile(delete=False)
    outfilename = outfile.name

    Entrez.email = "anatoly.urba[at]gmail.com"

    gis = []

    print("Scanning input.")
    with open(argparser.infile) as f:
        for line in f.readlines():
            if line.startswith("gi|"):
                gi, _t, descr = line.split('\t')
                if descr.strip() == "":
                    gi = gi.split('|')[1]
                    gis.append(gi)


    ann = {}
    gis = list(set(gis))
    print("%i ids to convert" % len(gis ))
    print("Fetching data.")

    for ggis in grouper(1024, gis):
        data = Entrez.esummary(db="protein", id = ",".join(ggis))
        entrez = Entrez.read(data)
        for e in entrez:
            ann[e['Gi']] = e['Title']
        print("%i ids fetched" % len(ann))

    print("Addinng annotation.")
    with open(argparser.infile) as f:
        for line in f.readlines():
            if line.startswith("gi|"):
                gi, _t, descr = line.split('\t')
                if descr.strip() == "":
                    gi_ = gi.split('|')[1]
                    descr = ann.get(int(gi_), "")
                line = '\t'.join([gi, _t, descr, '\n'])
            outfile.writelines(line)

    outfile.close()

    print("Done.")
    if argparser.nobk:
        os.unlink(argparser.infile)
    else:
        print("Backup file saved.")
        os.rename(argparser.infile, argparser.infile+".bk")

    shutil.move(outfilename, argparser.infile)


if __name__ == '__main__':
    main()