#!/usr/bin/python2
import sys
import os
import tempfile
from argparse import ArgumentParser

def main():
    argparser = ArgumentParser(description="Add gi annotation to peptblast output")
    argparser.add_argument("infile", type = str,  help = "Peptblast output file.")
    argparser = argparser.parse_args()

    outfile = tempfile.NamedTemporaryFile(delete=False)
    outfilename = outfile.name

    gis = []

    with open(argparser.infile) as f:
        for line in f.readlines():
            if line.startswith("gi|"):
                gi, _t, descr = line.split('\t')
                if descr.strip() == "":
                    gi = gi.split('|')[1]
                    gis.append(gi)

    print gis


    print outfilename





if __name__ == '__main__':
    main()