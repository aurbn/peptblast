#!/usr/bin/python2
from __future__ import print_function
import os
import sys
from subprocess import call
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline
from argparse import ArgumentParser

def simstr(alingment):
    res = []
    for i, j in zip (*alingment):
        if i == j:
            res.append("1")
        else:
            res.append("0")
    return "".join(res)

def encode(input_string):
    count = 1
    prev = ''
    lst = []
    for i,character in enumerate(input_string, 1):
        if character != prev:
            if prev:
                entry = (prev, count, ii)
                lst.append(entry)
                #print lst
            count = 1
            prev = character
            ii = i
        else:
            count += 1
    else:
        entry = (character,count, i)
        lst.append(entry)
    return lst

def thrids(z):
    return zip(z, z[1:], z[2:])

argparser = ArgumentParser(description="Peptide/protein blast helper tool")
argparser.add_argument('--db', type=str, required=True, help="Database.fasta")
argparser.add_argument('--pep', type=str, required=True, help="Peptides.fasta")
argparser.add_argument('--leftmin', type=int, required=False, default=1, help="Minimal length of left fragment")
argparser.add_argument('--rightmin', type=int, required=False, default=1, help="Minimal length of right fragment")
argparser.add_argument('--summin', type=int, required=False, default=6, help="Minimal total mathed length")
argparser.add_argument('--gapmax', type=int, required=False, default=3, help="Minimal gap size")
argparser.add_argument('-s', action='store_true', required=False, help='Use blast-short version.')
argparser.add_argument('-S', action='store_true', required=False, help='Suppress 5+ hits in 2 file.')
argparser = argparser.parse_args()

if not os.path.exists(argparser.db+".pin"):
    call("makeblastdb -dbtype prot -in "+argparser.db, shell=True)

if not os.path.exists(argparser.pep+".xml"):
	blast_cl = NcbiblastpCommandline(query = argparser.pep,
                                 db = argparser.db,
                                 task = "blastp-short" if argparser.s else "blastp",
                                 outfmt = 5,
                                 evalue = 0.1,
                                 num_threads = 10,
                                word_size = 2,
                                 out = argparser.pep+".xml"
                                )
	stdout, stderr = blast_cl()

if argparser.s:
    count = 0
    out = open(argparser.pep+'.1.txt', 'w')
    qresults = SearchIO.parse(argparser.pep+".xml", 'blast-xml')
    for qresult in qresults:
        for hit in qresult:
            for hsp in hit:
                if((hsp.aln_span == 5 and (hsp.gap_num == 0) and
                   (hsp.aln_span == hsp.ident_num)) or (hsp.aln_span > 5)):
                    count += 1
                    print(hsp, file = out)
                    print('', file = out)
    out = open(argparser.pep+'.2.txt', 'w')
    qresults = SearchIO.parse(argparser.pep+".xml", 'blast-xml')
    for qresult in qresults:
        for hit in qresult:
            for hsp in hit:
                if(hsp.aln_span >= 5 and (hsp.gap_num == 0) and
                   (hsp.aln_span == hsp.ident_num)):
                    count += 1
                    print(hsp, file = out)
                    print('', file = out)
else:
    count = 0
    out = open(argparser.pep+'.1.txt', 'w')
    qresults = SearchIO.parse(argparser.pep+".xml", 'blast-xml')
    for qresult in qresults:
        for hit in qresult:
            for hsp in hit:
                for v, c, p in encode(simstr(hsp.aln)):
                    if v=="1" and c>=5:
                        print( "5+ hit @ %i" % p, file = out)
                        print(hsp.aln[:,p-20:p+20], file = out)
                        print(simstr(hsp.aln)[p-20: p+20].replace("0", " ").replace("1", "+"),
                              file = out)
                        print('', file = out)
        

    out = open(argparser.pep+'.2.txt', 'w')
    qresults = SearchIO.parse(argparser.pep+".xml", 'blast-xml')
    for qresult in qresults:
        for hit in qresult:
            for hsp in hit:
                for t0, t1, t2 in thrids(encode(simstr(hsp.aln))):
                    if t0[0] == "1":
                        assert t2[0] == "1"
                        assert t1[0] == "0"
                        if t0[1] >= argparser.leftmin and t2[1] >= argparser.rightmin and\
                             (t0[1]+t2[1]) >= argparser.summin and t1[1] <= argparser.gapmax:
                            if not ( argparser.S and (t0[1] >= 5 or t2[1] >= 5 )):
                                p = t0[2]
                                print ("%i/%i/%i hit @ %i" % (t0[1], t1[1], t2[1], p), file = out)
                                print(hsp.aln[:,p-20:p+20], file = out)
                                print(simstr(hsp.aln)[p-20: p+20].replace("0", " ").replace("1", "+"),
                                      file = out)
                                print('', file = out)
