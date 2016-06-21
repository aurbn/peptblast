from __future__ import print_function
import os
import sys
from subprocess import call
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline
from argparse import ArgumentParser


argparser = ArgumentParser(description="Peptide/protein blast helper tool")
argparser.add_argument('--db', type=str, required=True, help='"Database.fasta"')
argparser.add_argument('--pep', type=str, required=True, help='"Peptides.fasta"')
argparser.add_argument('--s', action='store_true', required=False, help='Use blast-short version.')
argparser = argparser.parse_args()

if not os.path.exists(sys.argv[1]+".pin"):
    call("makeblastdb -dbtype prot -in "+argparser.db, shell=True)


blast_cl = NcbiblastpCommandline(query = argparser.pep,
                                 db = argparser.db,
                                 task = "blastp-short" if argparser.s else "blastp",
                                 outfmt = 5,
                                 evalue = 0.1,
                                 num_threads = 10,
                                word_size = 2,
                                 out = argparser.db+".xml"
                                )
stdout, stderr = blast_cl()

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
