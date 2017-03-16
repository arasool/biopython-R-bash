#!/usr/bin/env python
#creating a local blast for queries.
#blast_for_cgi.py script should be in the directory
#the blast script is not embedded here.
'''
Pythonic BLAST Web Interface
Write a program named runblast.cgi, which contains a web form that accepts the following
inputs:
1) Sequence - A FASTA sequence
2) Seqtype - Is the query sequence nt or protein
3) Blastprogram – which version of blast to use
4) Database – which database to search.
The output is the “hit id” of all sequences found.
This cgi script can be designed in two ways:
1) standalone cgi-script that executes your previous blast program OR
2) the blast program is embedded in the script.
There are advantages and disadvantages for both methods but feel free to pick the one you
think is easiest for you.
'''

print "Content-type: text/plain\r\n\r\n"


import cgi

# Create instance of FieldStorage
args = cgi.FieldStorage()

print "%s" % args

print "Sequence:", args.getvalue("sequence")
print "Seqtype:", args.getvalue("seqtype")
print "Blastprogram:", args.getvalue("blastprog")
print "Database:", args.getvalue("db")
print "E-value:", args.getvalue("evalue")


import blast

# write the input sequence to a file first
with open('tmp.fa', 'w') as f:
    f.write(args.getvalue("sequence"))
    f.write("\n")

# call run blast
ret, blastoutput = blast.run_blast('tmp.fa', args.getvalue("seqtype"))

print "blast return value:", ret
# print output
blast.show_alignment_scores(blastoutput, float(args.getvalue("evalue")))


#print "<html><body>Hello world from Python!</body></html>"
