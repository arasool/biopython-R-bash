###############################################################################
##this is question 1 on the midterm. 
##our input is a gff file, and our goal is to create a parser that would provide a list of annotated genes.
##I am going to import argparse, it is module makes it easy to write user-friendly command-line interfaces.
###############################################################################

'''In the first part of the midterm your task is to write a parser that will find and output a
list of genes that are annotated in a given region of the genome (on both strands).
There are cases when Biologists identify QTL (Quantitative Trait Loci) region which they
believe is responsible for a certain phenotype that they interested in. They often want
to identify all the genes that are present in that region so they can hypothesize which of
the genes is responsible for the phenotype. This script would come in handy so that they
can retrieve the genes that are present in that QTL region.
The file format that contains gene locus annotations is called a GFF file. There may be
modules written to parse them, so you are more than welcome to use them. However,
the format is a simple tab delimited format containing 9 columns, so it may easier to
simply write your own parser. The nine columns of a GFF file are :
1) Reference sequence – for example a chromosome name
2) Source – who created the annotations
3) Feature – genomic features such as gene, exons, CDS, etc. Notice they are
hierarchical. For our task we are only interested in the “gene” features.
4) Start – Start position of the feature on the Reference Sequence
5) End – Stop position of the feature on the Reference Sequence
6) Score – Such as Blast score, if there isn’t one you will see a “.”
7) Strand – Positive or Negative
8) Phase – Used for CDS to explain which of the three phases is being
translated.
9) Annotation – Details of the feature such as name and parent feature. For this
example we are only interested in “Name”
Write a Python script named midterm1.py that accepts five commandline options
“-i” => path to the GFF file
“-c” => Chromosome name
“-s” => Start coordinate of the region
“ -e ” => End coordinate of the region
“ -o ” => Where the output should be saved. If it is not provided, then it prints to
the screen.
'''


import argparse
##function parse_line  keep the file tab separated but get rid of whitespace character.
def parse_line(line):
    a = line.strip().split("\t")
##since we are interested in column 9, we are going to the name of the gene from there.
##k,v are key and value for the dictionary created here.
    annotation_pairs = a[8].split(";")
    annotation = {}
    for pair in annotation_pairs:
        if pair: # ignore empty
            k,v = pair.split("=")
            annotation[k] = v
##defining the header, columns in the file.
    obj = {
        'chromosome': a[0],
        'source': a[1],
        'feature': a[2],
        'start': int(a[3]),
        'end': int(a[4]),
        'score': a[5],
        'strand': a[6],
        'phase': a[7],
        'annotation': annotation,
    }
    return obj

##is_match function explains that we are only interested in 'gene' from the chromosome and we also want to 
##be able to specify start and end location on the chromosome in the command-line option.

def is_match(obj, args):
    conditions = [
        obj['feature'] == 'gene',
        obj['chromosome'] == args.chromosome,
        obj['start'] > args.start,
        obj['end'] < args.end
    ]
    return all(conditions)

###this module creation was new to me, and took some research into argparse  to declare all the command-line options.
##help=defines teh command line options.
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input', help='input GFF file')
    parser.add_argument('-c', dest='chromosome', help='chromosome name')
    parser.add_argument('-s', dest='start', type=int, help='start coordinate')
    parser.add_argument('-e', dest='end', type=int, help='end coordinate')
    parser.add_argument('-o', dest='output', help='output file')
    parser.add_argument('-f', dest='fasta', type=bool, help='True if output to fasta file')##for the extra credit, not sure if I will be able to get to it.
    parser.add_argument('-a', dest='fastafile', help='path to fasta reference file')
    args = parser.parse_args()
##the user has the option to either print the output or write it into an outfile.
    if args.output:
        outfile = open(args.output, 'w')

    for line in open(args.input):
        obj = parse_line(line)
        if is_match(obj, args):
            if not args.output:
                print obj['annotation']['Name']
            else:
                outfile.write("%s\n" % obj['annotation']['Name'])

    if args.output:
        outfile.flush()
        outfile.close()


