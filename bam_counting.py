## bam.counting.py

##commandline input arguments -g - .gff file, -b - .bam file, -d - database output filename
##sort and index bam using pysam

'''Write a program named bam counting.py which accepts input arguments (from the commandline)
of:
a GFF le (
ag \-g")
a BAM le (
ag \-b")
a sqlite3 database lename (
ag \-d", this is an output lename)
(use optparse/argparse to handle input options, make sure your code only runs when the code is
executed i.e. run bam counting.py or python bam counting.py) The program should:
1. load the contents of the input GFF le (-g)
2. sort and index the input BAM le using pysam
3. load the sorted and indexed BAM le
4. for every gene, extract the number of counts present in the BAM le (-b)
5. create an sqlite3 database (lename from -d) with a table named \Genes" with the schema
ID (integer, as the primary key)
GeneName (text)
seqname (text)
source (text)
start (integer)
end (integer)
score (text)
strand (text)
frame (text)
counts (integer)
6. add every gene found in the GFF le to the database
note: \ID" will be set by the database (or you can set it yourself), \counts" is the number
of counts found between the gene start and end indices, and the remaining information can be
extracted from the GFF le
'''


############################################
#!/usr/bin/python
import argparse
import sys
import pysam
import sqlite3
import os


##function parse_line parses the file, while keeping the file tab separated but get rid of whitespace character.

def parse_line(line):
    a = line.strip().split("\t")

##k,v are key and value for the dictionary created here.    

    annotation_pairs = a[8].split(";")
    annotation = {}
    for pair in annotation_pairs:
        if pair: # ignore empty
            k,v = pair.split("=")
            annotation[k] = v
##defining the header, columns in the file.
    obj = {
        'seqname': a[0],
        'source': a[1],
        'feature': a[2],
        'start': int(a[3]),
        'end': int(a[4]),
        'score': a[5],
        'strand': a[6],
        'frame': a[7],
        'annotation': annotation,
    }
    return obj

##Return a list of (counts, length) tuples for each gene feature.
##looking into bam file, (bam.references) I observed chr 1-5 then ChrC and ChrM.
##Return a list of (counts, length) tuples for each gene feature.
def cmap(seqname):
    if seqname == 'ChrC':
        return 'chloroplast'
    elif seqname == 'ChrM':
        return 'mitochondria'
    else:
        return seqname #chromosome number
##get_count is taken from the lab notes, just modified. It is taking the seqname (chloroplast or mitochonria and), start and end
##position of the chromosome and calculating total number of reads across chromosomes.

## "counts" is the number of counts found between the gene start and end indices.
def get_count(gff, bam):
    return len(list(bam.fetch(cmap(gff['seqname']), gff['start'], gff['end'])))

##copied the structure from sqlite tutorial online on c.execute. these ??? are like place holders, so 9 of those
##for 9 attributes we need to look into. this schema is given in the HW7 pdf.

def insert_row(c, gff, bam):
    c.execute('''INSERT INTO Genes (GeneName, seqname, source, start, end, score, strand, frame, counts)
        VALUES (?,?,?,?,?,?,?,?,?)''',
        (
            gff['annotation']['Name'],
            gff['seqname'],
            gff['source'],
            gff['start'],
            gff['end'],
            gff['score'],
            gff['strand'],
            gff['frame'],
            get_count(gff, bam)   #get_count is the function shown at line 53.
        )
    )

##using argparse to create the main method. Modified from midterm question 1.
##it has command-line option, with the type of file which would be accepted with that option.

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', dest='input', help='input GFF file')
    parser.add_argument('-b', dest= 'bam', help= 'bam file')
    parser.add_argument('-d', dest='output', help='output file')
    args = parser.parse_args()

    ##verifying the kind of argument that will be accepted.(A bam file, gff as input and the db file that is created as output of first two)

    if not args.bam or not args.input or not args.output:
        exit("Arguments missing.")
##got help on this part, basically using os module, determine the outfile and remove it for ongoing process until
## the whole process is completed.
    if os.path.exists(args.output):
        os.remove(args.output)

##using pysam to sort the bam file, taken this part from lab 9.
##Here we query a region and extract information (such as the counts). You can provide the reference name, start and end of the genomic region.
    pysam.sort(args.bam, 'sorted')
    pysam.index('sorted.bam')
    bam = pysam.AlignmentFile('sorted.bam', 'rb')

    conn = sqlite3.connect(args.output)

    #defining c
    c = conn.cursor()

    ##condition 5 on the homework pdf file.(Create an sqlite3 database (Filename from -d) with a table named \Genes" with the schema)

    c.execute('''CREATE TABLE Genes
        (id integer primary key,
        GeneName text,
        seqname text,
        source text,
        start integer,
        end integer,
        score text,
        strand text,
        frame text,
        counts integer)''')



    for line in open(args.input):
        gff = parse_line(line)
        if gff['feature'] == 'gene':
            insert_row(c, gff, bam)


    # Save (commit) the changes
    conn.commit()

    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()