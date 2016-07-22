#!/usr/bin/python2
from __future__ import print_function
import os
import sys
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
from argparse import ArgumentParser
import itertools
import shutil
import warnings
import StringIO
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO, SeqIO

TMP_DIR = "tmp_blast"


def getfastafiles(dir):
    return (os.path.join(dir, f) for f in os.listdir(dir)
            if os.path.isfile(os.path.join(dir, f)) and
            f.endswith(".fasta"))

def simstr(alingment):
    res = []
    for i, j in zip(*alingment):
        if i == j:
            res.append("1")
        else:
            res.append("0")
    return "".join(res)


def encode(input_string):
    count = 1
    prev = ''
    lst = []
    for i, character in enumerate(input_string, 1):
        if character != prev:
            if prev:
                entry = (prev, count, ii)
                lst.append(entry)
            count = 1
            prev = character
            ii = i
        else:
            count += 1
    else:
        entry = (character, count, i)
        lst.append(entry)
    return lst


def thrids(z):
    return zip(z, z[1:], z[2:])


def grouper(n, iterable):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk

def process_blast_output(file, simple):
        if simple:
            count = 0
            out1 = StringIO.StringIO()
            qresults = SearchIO.parse(file, 'blast-xml')
            for qresult in qresults:
                for hit in qresult:
                    for hsp in hit:
                        if ((hsp.aln_span == 5 and (hsp.gap_num == 0) and
                                 (hsp.aln_span == hsp.ident_num)) or (hsp.aln_span > 5)):
                            count += 1
                            print(hsp, file=out1)
                            print('', file=out1)
            out2 = StringIO.StringIO()
            qresults = SearchIO.parse(file, 'blast-xml')
            for qresult in qresults:
                for hit in qresult:
                    for hsp in hit:
                        if (hsp.aln_span >= 5 and (hsp.gap_num == 0) and
                                (hsp.aln_span == hsp.ident_num)):
                            count += 1
                            print(hsp, file=out2)
                            print('', file=out2)
    else:
        count = 0
        out1 = StringIO.StringIO()
        qresults = SearchIO.parse(file, 'blast-xml')
        for qresult in qresults:
            for hit in qresult:
                for hsp in hit:
                    for v, c, p in encode(simstr(hsp.aln)):
                        if v == "1" and c >= 5:
                            print("5+ hit @ %i" % p, file=out1)
                            print(hsp.aln[:, p - 20:p + 20], file=out1)
                            infs = simstr(hsp.aln)
                            infs = infs[:p - 1] + "#" * c + infs[p + c + 22:]
                            infs = infs[p - 20: p + 20].replace("0", " ").replace("1", "+")
                            print(infs, file=out1)
                            print('', file=out1)

        out2 = StringIO.StringIO()
        qresults = SearchIO.parse(argparser.pep + ".xml", 'blast-xml')
        for qresult in qresults:
            for hit in qresult:
                for hsp in hit:
                    for t0, t1, t2 in thrids(encode(simstr(hsp.aln))):
                        if t0[0] == "1":
                            assert t2[0] == "1"
                            assert t1[0] == "0"
                            if t0[1] >= argparser.leftmin and t2[1] >= argparser.rightmin and \
                                            (t0[1] + t2[1]) >= argparser.summin and t1[1] <= argparser.gapmax:
                                if not ( argparser.S and (t0[1] >= 5 or t2[1] >= 5 )):
                                    p = t0[2]
                                    print("%i/%i/%i hit @ %i" % (t0[1], t1[1], t2[1], p), file=out2)
                                    print(hsp.aln[:, p - 20:p + 20], file=out2)
                                    print(simstr(hsp.aln)[p - 20: p + 20].replace("0", " ").replace("1", "+"),
                                          file=out2)
                                    print('', file=out2)



def main():
    argparser = ArgumentParser(description="Peptide/protein blast helper tool")
    argparser.add_argument('--db', type=str, required=True, help="Database.fasta")
    argparser.add_argument('--pep', type=str, required=True, help="Peptides.fasta")
    argparser.add_argument('--leftmin', type=int, required=False, default=1, help="Minimal length of left fragment")
    argparser.add_argument('--rightmin', type=int, required=False, default=1, help="Minimal length of right fragment")
    argparser.add_argument('--summin', type=int, required=False, default=6, help="Minimal total mathed length")
    argparser.add_argument('--gapmax', type=int, required=False, default=3, help="Minimal gap size")
    argparser.add_argument('--threads', type=int, required=False, default=4, help="Blast threads")
    argparser.add_argument('-s', action='store_true', required=False, help='Use blast-short version.')
    argparser.add_argument('-S', action='store_true', required=False, help='Suppress 5+ hits in 2 file.')
    argparser = argparser.parse_args()

    if not os.path.exists(argparser.db + ".pin"):
        subprocess.call("makeblastdb -dbtype prot -in " + argparser.db, shell=True)

    if os.path.exists(TMP_DIR):
        shutil.rmtree(TMP_DIR)

    os.mkdir(TMP_DIR)
    # Determine number of records in fasta file
    proc = subprocess.Popen("grep '>' " + argparser.pep + " | wc -l", shell=True,
                            stdout=subprocess.PIPE).communicate()[0]
    proc = int(proc)
    assert proc > 0 and proc < 100

    chunksize = proc / argparser.threads
    inp = SeqIO.parse(open(argparser.pep), "fasta")
    # Write equal number of records for each thread
    for i, chunk in enumerate(grouper(chunksize, inp)):
        SeqIO.write(chunk, open(os.path.join(TMP_DIR, str(i)+".fasta"), "w"), "fasta")

    # Run `threads` number of blast commands in parallel
    processes = set()
    for file in getfastafiles(TMP_DIR):
        blast_cl = NcbiblastpCommandline(query=file,
                                         db=argparser.db,
                                         task="blastp-short" if argparser.s else "blastp",
                                         outfmt=5,
                                         evalue=100,
                                         num_threads=1,
                                         word_size=2,
                                         out= file.replace(".fasta", ".xml")
        )
        print(str(blast_cl))
        processes.add(subprocess.Popen(str(blast_cl), shell=True))
        if len(processes) >= argparser.threads:
            os.wait()
            processes.difference_update([
                p for p in processes if p.poll() is not None])



if __name__ == "__main__":
    main()