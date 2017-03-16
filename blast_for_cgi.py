
import os.path
import subprocess
import sys

from Bio.Blast import NCBIXML


def show_alignment_scores(blast_output, evalue=0.01):
    # blast output is an XML
    records = NCBIXML.parse(open("/tmp/" + blast_output))

    #  Print header
    # print "QueryID\tSbjctAcc\tQueryCover(%)\tTotalScore(bits)"

    for record in records:
        query_len = float(record.query_length)
        query_id = record.query_id
        for alignment in record.alignments:
            print alignment.hit_id
            for hsp in alignment.hsps:
                if hsp.expect <= evalue:
                    hsp_flag = True
                    total_score += hsp.bits
                    query_cover = query_cover.union(range(hsp.query_start,hsp.query_end+1))
            
             Only print line if there is at least one significant HSP in alignment
            if hsp_flag:
                print '\t'.join(map(str, [query_id, alignment.accession, round((len(query_cover)/query_len)*100,2), round(total_score,2)]))


def run_blast(input_file, sequence_type='prot', prog='blastp'):

    dbprefix = '/usr/local/share/Blast/db/'  #within cygnus

    # figure out xml filename
    filename = input_file.rsplit(".", 1)[0] # get filename without extension
    output_file = "%s.xml" % filename

    # NOTE: would have to download all nr.##.tar.gz and nt.##.tar.gz at some NIH site
    #       for this to work. It's >20 files that are all 500MB-1GB large, in case we want to run this locally without cygnus
    # create the command
    if sequence_type == 'prot':
        cmd = '%s -query %s -db %snr -outfmt 5 -out /tmp/%s' % (prog, input_file, dbprefix, output_file)
    elif sequence_type == 'na':
        cmd = '%s -query %s -db %snt -outfmt 5 -out /tmp/%s' % (prog, input_file, dbprefix, output_file)
    print "Running:", cmd
    return subprocess.call(cmd.split(), stderr=subprocess.STDOUT), output_file



def main():
    # set defaults
    sequence_type = 'prot'
    evalue = 0.01

    # get args
    input_file = sys.argv[1]
    if len(sys.argv) == 4:
        # pythonfile inputfile na/prot evalue
        sequence_type = sys.argv[2]
        evalue = float(sys.argv[3])
    elif len(sys.argv) == 3:
        # pythonfile inputfile (evalue OR sequence type)
        if sys.argv[2] in ('prot', 'na'):
            sequence_type = sys.argv[2]
        else:
            evalue = float(sys.argv[2])

    retvalue, output_file = run_blast(input_file, sequence_type)
    show_alignment_scores(output_file, evalue)




if __name__ == '__main__':
    main()


# $ python blast.py PCYB DBP.fa prot 1
# $ python blast.py PCYB DBP.fa prot
# $ python blast.py PCYB DBP.fa 1e-5
# $ python blast.py PCYB DBP.fa
# $ python blast.py RBP1a.fa na
# $ python blast.py RBP1a.fa na 0.001
