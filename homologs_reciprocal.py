############################################################
##This is Question 2 on midterm. The purpose of this script is to find 
##identify all reciprocal best matches from the two different species.
##I need blast tool from Biopython, and I will be using XML files for input containing Blast outputs.
##For this question, we are only interested in maximum bit score, and the highest bit score gets replaced
## each time there is a higher value found during search in each alignment.
##just like the Blast tool site, we have query sequence that we provide and alignment is provided against it.
##so to find the reciprocal best match, one species' genomic information becomes query and it aligned against second specie.
##then second specie becomes query and it is aligned against the first.
##variable hitdef is used in place of HSP: high scoring pair.
##The output may provide the max bit score and the alignment pair that got the score.
#################################################################
'''Identifying Homologs based on reciprocal best BLAST hit.
Write a Python script and name it midterm2.py that takes as input two different Blast
outputs in xml format and identifies all reciprocal best matches from the two different
species. For this script, simply use the total bit-score to score your matches. If multiple
matches have the same score, then they will all be considered best hits (multiple
homologs).

'''



import sys
from Bio.Blast import NCBIXML


def is_reciprocal(query, hitdef, B2):
    for br2 in B2:
        if br2.query == hitdef:
            maxbitscore, maxhitdefs = compute_max_bitscore(br2)
            #print "Q:", query, "HD:", maxhitdefs
            return query in maxhitdefs

##this function defines how to find the max bit score and keep recording the maximum score that is found.
def compute_max_bitscore(br):
    maxbitscore = 0
    maxhitdefs = []

    for alignment in br.alignments:
        bitscore = sum(hsp.bits for hsp in alignment.hsps)
        if bitscore > maxbitscore:
            maxbitscore = bitscore
            maxhitdefs = [alignment.hit_def]
        elif bitscore == maxbitscore:
            maxhitdefs.append(alignment.hit_def)

    return maxbitscore, maxhitdefs

##homologs.txt is declared as the output file. All the results will be printed in this file.
def compare(B1, B2):
    with open('homologs.txt', 'w') as f:

        for br1 in B1:
            maxbitscore, maxhitdefs = compute_max_bitscore(br1)
            # print maxbitscore, maxhitdefs
##the result will be written in tab separated format.
            for hitdef in maxhitdefs:
                if is_reciprocal(br1.query, hitdef, B2):
                    line = "%s\t%s\n" % (br1.query, hitdef)
                    f.write(line)
                    print line,

## need two arguments, two alignment files to find reciprocal best hit. 
if __name__ == '__main__':
    if len(sys.argv) < 3:
        exit("2 arguments required.")

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    # print file1, file2

    B1 = NCBIXML.parse(open(file1))
    B2 = NCBIXML.parse(open(file2))
##calling out the functions defined above.
    compare(B1, B2)
