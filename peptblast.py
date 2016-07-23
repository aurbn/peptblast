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
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO, SeqIO

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
                    if ((hsp.aln_span == 5 and (hsp.gap_num == 0) and
                             (hsp.aln_span == hsp.ident_num)) or (hsp.aln_span > 5)):
                        yield (str(hsp), None)
                        yield ('\n\n', None)

                for hsp in hit:
                    if (hsp.aln_span >= 5 and (hsp.gap_num == 0) and
                            (hsp.aln_span == hsp.ident_num)):
                        yield (None, str(hsp))
                        yield (None, '\n\n')
    else:
        for qresult in qresults:
            for hit in qresult:
                for hsp in hit:
                    for v, c, p in encode(simstr(hsp.aln)):
                        if v == "1" and c >= 5:
                            yield (format_alignment(hsp, p, c), None)
                for hsp in hit:
                    for t0, t1, t2 in thrids(encode(simstr(hsp.aln))):
                        if t0[0] == "1":
                            assert (t0[2] < t1[2] < t2[2])
                            assert t2[0] == "1"
                            assert t1[0] == "0"
                            if t0[1] >= argparser.leftmin and t2[1] >= argparser.rightmin and \
                                            (t0[1] + t2[1]) >= argparser.summin and t1[1] <= argparser.gapmax:
                                if not ( argparser.S and (t0[1] >= 5 or t2[1] >= 5 )):
                                    yield (None, format_alignment(hsp, t0[2], t0[1], t2[2], t2[1]))


def format_alignment(hsp, pos1, len1, pos2=None, len2=None, rightgap=10, maxlen=50):
    assert not (pos2 is None) ^ (len2 is None)
    qstr = str(hsp.query.seq.lower())
    hstr = str(hsp.hit.seq.lower())
    sstr = str(hsp.aln_annotation['similarity'].lower())
    lines = []

    if pos2 is None:  # 5+hit
        lines.append("%i length hit at %i\n" % (len1, pos1))
        l1 = len(qstr)
        qstr = qstr[:pos1-1]+qstr[pos1-1:pos1+len1-1].upper()+qstr[pos1+len1-1:]
        hstr = hstr[:pos1-1]+hstr[pos1-1:pos1+len1-1].upper()+hstr[pos1+len1-1:]
        sstr = sstr[:pos1-1]+sstr[pos1-1:pos1+len1-1].upper()+sstr[pos1+len1-1:]
        assert l1 == len(qstr)
    else: # complex hit
        lines.append("%i/%i/%i hit at %i\n" % (len1, pos2-(pos1+len1), len2, pos1))
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

    qstr = qstr + "\t" + hsp.query_id + "\t" + hsp.query_description
    hstr = hstr + "\t" + hsp.hit_id + "\t" + hsp.hit_description
    qstr += '\n'
    lines.append(qstr)
    sstr += '\n'
    lines.append(sstr)
    hstr += '\n\n'
    lines.append(hstr)

    return lines



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

    if os.path.exists(TMP_DIR):
        shutil.rmtree(TMP_DIR)
    os.mkdir(TMP_DIR)

    db_name = os.path.join(TMP_DIR, argparser.pep)

    subprocess.call("makeblastdb -dbtype prot -in " + argparser.db
                    + " -out " + db_name, shell=True)


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
    for file in getfiles(TMP_DIR, "fasta"):
        blast_cl = NcbiblastpCommandline(query=file,
                                         db=db_name,
                                         task="blastp-short" if argparser.s else "blastp",
                                         outfmt=5,
                                         evalue=100,
                                         num_threads=1,
                                         word_size=2,
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


    out1 = open(argparser.pep + "1.txt", "w")
    out2 = open(argparser.pep + "2.txt", "w")
    for file in getfiles(TMP_DIR, "xml"):
        print ("Reading " + file)
        for s1, s2 in process_blast_output(file, argparser.s, argparser):
            if s1 is not None:
                out1.writelines(s1)
            if s2 is not None:
                out2.writelines(s2)
    out1.close()
    out2.close()


if __name__ == "__main__":
    main()