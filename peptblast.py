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
import operator
import multiprocessing
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO, SeqIO

__version__ = "0.2"

TMP_DIR = "tmp_blast"


def getfiles(dir, extension):
    return (os.path.join(dir, f) for f in os.listdir(dir)
            if os.path.isfile(os.path.join(dir, f)) and
            f.endswith("."+extension))


def simstr(alignment):
    res = []
    for i, j in zip(*alignment):
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


def process_blast_output(file, simple, argparser):
    qresults = SearchIO.parse(file, 'blast-xml')
    if simple:
        for qresult in qresults:
            for hit in qresult:
                for hsp in hit:
                    if ((hsp.aln_span == argparser.cont and (hsp.gap_num == 0) and
                             (hsp.aln_span == hsp.ident_num)) or (hsp.aln_span > argparser.cont)):
                        yield ([str(hsp), "\n\n"], None, hsp.aln_span)

                for hsp in hit:
                    if (hsp.aln_span >= argparser.cont and (hsp.gap_num == 0) and
                            (hsp.aln_span == hsp.ident_num)):
                        yield (None, [str(hsp), "\n\n"], hsp.aln_span)
    else:
        for qresult in qresults:
            for hit in qresult:
                for hsp in hit:
                    for v, c, p in encode(simstr(hsp.aln)):
                        if v == "1" and c >= argparser.cont:
                            yield (format_alignment(hsp, p, c), None, c)
                for hsp in hit:
                    for t0, t1, t2 in thrids(encode(simstr(hsp.aln))):
                        if t0[0] == "1":
                            assert (t0[2] < t1[2] < t2[2])
                            assert t2[0] == "1"
                            assert t1[0] == "0"
                            if t0[1] >= argparser.leftmin and t2[1] >= argparser.rightmin and \
                               (t0[1] + t2[1]) >= argparser.summin and \
                                t1[1] <= argparser.gapmax:
                                if not (argparser.S and
                                        (t0[1] >= argparser.cont or t2[1] >= argparser.cont)):
                                    yield (None, format_alignment(hsp, t0[2], t0[1], t2[2], t2[1]),
                                           t0[1]+t2[1]-t1[1])


def format_alignment(hsp, pos1, len1, pos2=None, len2=None, rightgap=10, maxlen=100):
    assert not (pos2 is None) ^ (len2 is None)
    qstr = str(hsp.query.seq.lower())
    hstr = str(hsp.hit.seq.lower())
    sstr = str(hsp.aln_annotation['similarity'].lower())
    lines = []

    if pos2 is None:  # 5+hit
        lines.append("%i length hit at [%i:%i]\n" % (len1, pos1+hsp.query_start,
                                                     len1+pos1+hsp.query_start ))
        l1 = len(qstr)
        qstr = qstr[:pos1-1]+qstr[pos1-1:pos1+len1-1].upper()+qstr[pos1+len1-1:]
        hstr = hstr[:pos1-1]+hstr[pos1-1:pos1+len1-1].upper()+hstr[pos1+len1-1:]
        sstr = sstr[:pos1-1]+sstr[pos1-1:pos1+len1-1].upper()+sstr[pos1+len1-1:]
        assert l1 == len(qstr)
    else: # complex hit
        lines.append("%i/%i/%i hit at [%i:%i]\n" % (len1, pos2-(pos1+len1), len2, pos1+hsp.query_start,
                                                    pos2+len2+hsp.query_start))
        l1 = len(qstr)
        # TODO: fix bug with end-line alignments
        qstr = qstr[:pos1-1]+qstr[pos1-1:pos1+len1-1].upper()+qstr[pos1+len1-1:pos2-1]+qstr[pos2-1:pos2+len2-1].upper()+qstr[pos2+len2-1:]
        hstr = hstr[:pos1-1]+hstr[pos1-1:pos1+len1-1].upper()+hstr[pos1+len1-1:pos2-1]+hstr[pos2-1:pos2+len2-1].upper()+hstr[pos2+len2-1:]
        sstr = sstr[:pos1-1]+sstr[pos1-1:pos1+len1-1].upper()+sstr[pos1+len1-1:pos2-1]+sstr[pos2-1:pos2+len2-1].upper()+sstr[pos2+len2-1:]
        assert l1 == len(qstr)

    if pos1 > rightgap:
        qstr = qstr[pos1-rightgap:]
        hstr = hstr[pos1-rightgap:]
        sstr = sstr[pos1-rightgap:]
    if len(sstr) > maxlen:
        # TODO: check if match more than maxlen
        qstr = qstr[:maxlen]
        hstr = hstr[:maxlen]
        sstr = sstr[:maxlen]

    qstr = hsp.query_id + "\t[%i : %i]\t" % (hsp.query_start, hsp.query_end) + hsp.query_description + '\n' \
           + qstr
    hstr = hstr + \
           '\n' + hsp.hit_id + "\t[%i : %i]\t" % (hsp.hit_start, hsp.hit_end) + hsp.hit_description
    qstr += '\n'
    lines.append(qstr)
    sstr += '\n'
    lines.append(sstr)
    hstr += '\n\n'
    lines.append(hstr)

    return lines


def makeblastdb(infile, name):
     subprocess.call("makeblastdb -dbtype prot -in " + infile
                    + " -out " + name, shell=True)


def parallel_blast(argparser, db_name):
    # Determine number of records in fasta file
    nrec = subprocess.Popen("grep '>' " + argparser.pep + " | wc -l", shell=True,
                            stdout=subprocess.PIPE).communicate()[0]
    nrec = int(nrec)
    assert nrec > 0

    chunksize = nrec / argparser.threads
    chunksize = chunksize if chunksize else nrec
    inp = SeqIO.parse(open(argparser.pep), "fasta")

    # Write equal number of records for each thread
    for i, chunk in enumerate(grouper(chunksize, inp)):
        SeqIO.write(chunk, open(os.path.join(TMP_DIR, str(i)+".fasta"), "w"), "fasta")

    # Estimate required e-value
    if argparser.eval is None:
        filesize = os.path.getsize(argparser.db)
        e_est = filesize*1./(20**argparser.cont)
        print ("Estimated required evalue: %.2f" % e_est)
    else:
        e_est = argparser.eval

    # Run `threads` number of blast commands in parallel
    processes = set()
    for file in getfiles(TMP_DIR, "fasta"):
        blast_cl = NcbiblastpCommandline(query=file,
                                         db=db_name,
                                         task="blastp-short" if argparser.s else "blastp",
                                         outfmt=5,
                                         evalue=e_est,
                                         num_threads=1,
                                         word_size=argparser.wordsize,
                                         out= file.replace(".fasta", ".xml")
        )

        print ("Blasting " + file)
        processes.add(subprocess.Popen(str(blast_cl), shell=True))
        if len(processes) >= argparser.threads:
            os.wait()
            processes.difference_update([
                p for p in processes if p.poll() is not None])
    while processes:
        os.wait()
        processes.difference_update([
            p for p in processes if p.poll() is not None])


def read_blast_xml((file, argparser)):
    print("Reading " + file)
    results1 = []
    results2 = []
    for s1, s2, i in process_blast_output(file, argparser.printold, argparser):
        if s1 is not None:
            results1.append((s1, i))
        if s2 is not None:
            results2.append((s2, i))
    return (results1, results2)


def main():
    argparser = ArgumentParser(description="Peptide/protein blast helper tool")
    argparser.add_argument('--db', type=str, required=False, help="Database.fasta")
    argparser.add_argument('--pep', type=str, required=False, help="Peptides.fasta")
    argparser.add_argument('--cont', type=int, required=False, default=5, help="Minimal length of continuous fragment")
    argparser.add_argument('--leftmin', type=int, required=False, default=3, help="Minimal length of left fragment")
    argparser.add_argument('--rightmin', type=int, required=False, default=3, help="Minimal length of right fragment")
    argparser.add_argument('--summin', type=int, required=False, default=8, help="Minimal total mathed length")
    argparser.add_argument('--gapmax', type=int, required=False, default=3, help="Minimal gap size")
    argparser.add_argument('--threads', type=int, required=False, default=4, help="Blast threads")
    argparser.add_argument('--eval', type=float, required=False, help="Required e-value for blast, otherwise estimated")
    argparser.add_argument('--wordsize', type=int, required=False, default=2, help="Blast word size")
    argparser.add_argument('-s', action='store_true', required=False, help='Use blast-short instead of blast.')
    argparser.add_argument('--printold', action='store_true', required=False, help='Use old filtering/printing method.')
    argparser.add_argument('-S', action='store_true', required=False, help='Suppress 5+ hits in 2 file.')
    argparser.add_argument('--sort', action='store_true', required=False, help='Sort output.')
    argparser.add_argument('--keeptmp', action='store_true', required=False, help='Keep temporary files.')
    argparser.add_argument('--usetmp', type=str, required=False, default= None, help="Use existing temporary directory"
                                                                                                     "Do not perform blast.")

    argparser = argparser.parse_args()

    if argparser.usetmp is not None:
        if os.path.isdir(argparser.usetmp):
            argparser.pep = "pep"
            argparser.db = [f.replace(".pin", "") for f in os.listdir(argparser.usetmp) if ".pin" in f][0]
        else:
            print("Wrong tmp directory!")
            sys.exit(1)
    else:
        if (argparser.pep is None) or (argparser.pep is None):
            print ("Provide query and database or use tmp dir!")
            sys.exit(1)


    if os.path.exists(TMP_DIR) and not argparser.usetmp:
        shutil.rmtree(TMP_DIR)
    elif not argparser.usetmp:
        os.mkdir(TMP_DIR)

    if not argparser.usetmp:
        db_name = os.path.join(TMP_DIR, argparser.db)
        makeblastdb(argparser.db, db_name)
        parallel_blast(argparser, db_name)


    with open(argparser.pep + "_" + argparser.db + "_1.txt", "w") as out1, \
         open(argparser.pep + "_" + argparser.db + "_2.txt", "w") as out2:
        pool = multiprocessing.Pool(processes=argparser.threads)
        results1, results2 = zip(*pool.map(read_blast_xml, zip(getfiles(TMP_DIR, "xml"),
                                               itertools.repeat(argparser))))
        results1 = reduce(operator.add, results1, [])
        results2 = reduce(operator.add, results2, [])

        if argparser.sort:
            print ("Sorting results...")
            results1.sort(key = operator.itemgetter(1), reverse=True)
            results2.sort(key = operator.itemgetter(1), reverse=True)
        for r, i in results1:
            out1.writelines(r)
        for r, i in results2:
            out2.writelines(r)


    if not (argparser.keeptmp or argparser.usetmp):
        shutil.rmtree(TMP_DIR)

if __name__ == "__main__":
    main()