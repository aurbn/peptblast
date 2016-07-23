__author__ = 'urban'
from Bio import SearchIO

qrs = SearchIO.parse("tmp_blast/3.xml", "blast-xml")
for qr in qrs:
    for hit in qr:
        for hsp in hit:
            print(hsp)