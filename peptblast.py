from __future__ import print_function
import os
import sys
from subprocess import call
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline


if len(sys.argv) < 3:
    print("Usge: pepblast.py Database.fasta Peptides.fasta")
    sys.exit(1)

if not os.path.exists(sys.argv[1]+".pin"):
    call("makeblastdb -dbtype prot -in "+sys.argv[1], shell=True)


blast_cl = NcbiblastpCommandline(query = sys.argv[2],
                                 db = sys.argv[1],
                                 task = "blastp-short",
                                 outfmt = 5,
                                 evalue = 0.1,
                                 num_threads = 10,
                                word_size = 2,
                                 out = sys.argv[2]+".xml"
                                )
stdout, stderr = blast_cl()

count = 0
out = open(sys.argv[2]+'.1.txt', 'w')
qresults = SearchIO.parse(sys.argv[2]+".xml", 'blast-xml')
for qresult in qresults:
    for hit in qresult:
        for hsp in hit:
            if((hsp.aln_span == 5 and (hsp.gap_num == 0) and
               (hsp.aln_span == hsp.ident_num)) or (hsp.aln_span > 5)):
                count += 1
                print(hsp, file = out)
                print('', file = out)
out = open(sys.argv[2]+'.2.txt', 'w')
qresults = SearchIO.parse(sys.argv[2]+".xml", 'blast-xml')
for qresult in qresults:
    for hit in qresult:
        for hsp in hit:
            if(hsp.aln_span >= 5 and (hsp.gap_num == 0) and
               (hsp.aln_span == hsp.ident_num)):
                count += 1
                print(hsp, file = out)
                print('', file = out)
