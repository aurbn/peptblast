# peptblast #
A set of scripts providing convenient interface for finding peptide matches in proteins and/or peptide matches between non-orthologous proteins.

## peptblast.py ##
Script performs blast and generates two files: PEP_DB_1.txt contains continuous matches and PEP_DB_2.txt contains matches with a gap. 

### Parameters ###
    -h, --help           show help message and exit.  
    --db DB              Database in fasta format (required).  
    --pep PEP            Peptides or proteins file in fasta format (required).  
    --cont CONT          Minimal length of continuous fragment, default = 5.  
    --leftmin LEFTMIN    Minimal length of left fragment, default = 3.  
    --rightmin RIGHTMIN  Minimal length of right fragment, default = 3.  
    --summin SUMMIN      Minimal total matched length for fragmented hit, default = 8.  
    --gapmax GAPMAX      Minimal gap size, maximal gap size, default = 3.  
    --threads THREADS    Blast threads, default = 4.  
    --eval EVAL          Required e-value for blast, otherwise estimated as `approx.db_size/(20^cont)`.  
    --wordsize WORDSIZE  Blast word size, default = 2.  
    -s                   Use blast-short tool instead of blast.  
    --printold           Use old filtering/printing method.  
    -S                   Suppress continuous hots in _2.txt file.  
    --sort               Sort output, metric = length for contunous hits or `sum(matched)-length(gap)` for fragmented hits.  
    --keeptmp            Keep temporary files.  
    --usetmp USETMP      Use existing temporary directory for creating output. Do not perform blast.  

### Example command line ###
./peptblast.py --pep Targets.fasta --db UniProtKB.fasta --keeptmp --threads 24 --sort --cont 6

### Output ###
Each output hit stored as a 6-line text record with empty lines between them.
Example:

    6 length hit at [43:49]
    WP_002955589.1  [4 : 220]       MULTISPECIES: ABC transporter ATP-binding protein [Helicobacter]
    vfitgasgsGKSTILssffghlgvksghlnvfgvsmqkaskarinhlrrnigivfqdyklikewniernvmlpmvingykkevcksqvekllvhiklshk
    + + g sg+GKSTIL   f    + sg + + g  + + ++a    lr +ig+v qd  l  +  i  n+    v  g        +ve       +
    lalvgpsgaGKSTILrllfrfydissgciridgqdisqvtqa---slrshigvvpqdtvlfnd-tiadnirygrvtag------ndeveaaaqaagihda
    tr|A0A024R436|A0A024R436_HUMAN  [183 : 398]     ATP-binding cassette, sub-family B (MDR/TAP), member 6, isoform CRA_f OS=Homo sapiens GN=ABCB6 PE=3 SV=1

Lines description:

1. Contains 3 tab separated fields.  
 1.1. Length of continuous hit or `a/b/c` where a is left fragment length, c is right fragment length and b is gap length in fragmented hit.  
 1.2. Position of the hit in the blasts HSP.  
 1.3. Number of occurrences of db hit in the target (added by another script).  
2. Also contains 3 tab separated fields.  
 2.1. ID in the pep file.  
 2.2. Position of the hit in the pep file record.  
 2.3. Description of the hit in the pep file (optional).  
3. Hit sequence, hit is in UPPERCASE.  
4. Blasts midline, hit is in UPPERCASE.  
5. Hit sequence in db, hit is in UPPERCASE.  
6. As field 2 contains 3 tab separated fields.  
 6.1. ID in the db file  
 6.2. Position of the hit in the db file record.  
 6.3. Description of the hit in the pep file (optional).  

## merger_filer.py ##
Script performs merging _1.txt and _2.txt files produced by `peptblast.py` and basic filtration.

### Parameters ###

    -h, --help       show help message and exit.
    --1 ONE          continuous (_1.txt) hits file (required).
    --2 TWO          gaped (_2.txt) hits file (required).
    --filter FILTER  Only use less or equal cont occurrences (Can be applied only to `counter_filter.py` processed files).
    --best           Select only best match with query/hit combination based on metric (length for continuous hits a+c-b for gaped hits), default = False.
    --out OUT        Output file, default = merged.txt.

### Example command line ###
./merger_filter.py --1 PEP_DB_1.txt --2 PEP_DB_2.txt --best --out merged.txt

## gi_annotator.py ##
Adds GI annotations (fields 3) in lines 2 or 6 if not present. Modifies file in-place.

### Parameters ###

    infile      Peptblast output file (required).
    -h, --help  show help message and exit.
    --nobk      Do not save backup file.

### Example command line ###
./gi_annotator.py merged.txt

## counter_filter.py ##
Script counts queries and hits occurrences and performs filtration. Modifies file in-place.

### Parameters ###

    infile                Peptblast output file (required).
    -h, --help            show help message and exit.
    --count {query,hit}   What id to count: `query` or `hit`.
    --nosort              Don't sort and modify input file by occurencies.
    --table TABLE         Create output table of occurrences.
    --nobk                Do not save backup file.
    --taxfilter TAXFILTER Allow only hits with these taxons separated by commas (if organism unidentified the record is leaved).
    --organisms ORGANISMS   Create output table for organism frequencies.
    --le LE               Print only hits occurred less or equal times

### Example command line ###
./counter_filter.py merged.txt --count query --table queries.txt --taxfilter bacteria --organisms organinsms.txt --le 10

## descriptionfilter.py ##
Filter entries by description of query and hit. Delete records with fraction of common words in description >'minfcommon'.

### Parameters ###

    infile                  Peptblast output file.
    -h, --help              show this help message and exit
    --wordlen WORDLEN       Minimal length of a word to consider
    --minfcommon MINFCOMMON Minimal proportion of common words to all words in descriptions to filter out a record
    --out OUT               Output file.

### Example command line ###
./descriptionfilter.py --wordlen 3 --minfcommon 0.5

## freqsplitter.py ##
Splits input file according to query frequencies: `head` least frequent entries into infile_head, `tail` most frequent entries into infile_tail and rest into infile_mid.

### Parameters ###
    infile      Peptblast output file.
    head        Size of 'head'
    tail        Size of 'tail'
    
### Example command line ###
./freqsplitter.py merged.txt 50 50

# Typical workflow #
peptblast.py -> merger_filter.py -> gi_annotator.py -> descriptionfilter.py -> counter_filter.py -> freqsplitter.py


