#!/usr/bin/env python3

import subprocess
import io
import argparse


# master function. Takes in command line arguments and calls all necessary functions for the program to run
def annotate_genome(fasta_input, confidence_interval, seqname):
    remove_file("blast_fasta.fasta")
    remove_file("revcomp_fasta.fasta")

    # calling confidence_generation function and saving output to variable count. represents the number of codons needed without a stop codon to constitute an ORF
    count = confidence_generation(confidence_interval)
    gff_lines = []
    print("beginning searching for ORFs in reading frames")
    try:
        with open(fasta_input, "r") as fasta_handle:
            gene_search(fasta_handle, count, seqname, 0, gff_lines)
            gene_search(fasta_handle, count, seqname, 1, gff_lines)
            gene_search(fasta_handle, count, seqname, 2, gff_lines)
            line_lengths = write_revcomp(fasta_handle)
        with open("revcomp_fasta.fasta", "r") as fasta_handle:
            gene_search(fasta_handle, count, seqname, 3, gff_lines, line_lengths)
            gene_search(fasta_handle, count, seqname, 4, gff_lines, line_lengths)
            gene_search(fasta_handle, count, seqname, 5, gff_lines, line_lengths)
    except IOError as error:
        print("error reading input fasta", IOError)
    print("finished searching for ORFs in reading frames")
    blast()
    parse_blast_output(gff_lines)
    gff_writer(gff_lines)


# takes in command line argument confidence value and returns the value of the number of codons needed before hitting a stop codon in order to consider a string of nucleotides an ORF
def confidence_generation(confidence_interval):
    probability = 61/64
    confidence_interval = confidence_interval
    continuous = probability
    count = 1
    while(1-continuous < confidence_interval):
        count += 1
        continuous = continuous * probability
    return count

# takes the opened fasta input file and reverses and complements the entire thing and writes this to a new file called revcomp_fasta.fasta
def write_revcomp(fasta_handle):
    fasta_handle.seek(0)
    line_string = ""
    line_lengths = []
    prev_line_length = 0
    for line in fasta_handle:
        # if it is a header line we will write the previous line but revcomped, and then the header as normal
        if line[0] == ">":
            try:
                with open("revcomp_fasta.fasta", "a") as revcomp_handle:
                    if line_string != "":
                        line_lengths.append(len(line_string)+prev_line_length)
                        prev_line_length = len(line_string)
                        revcomp_handle.write(line_string[::-1])
                        line_string = ""
                    revcomp_handle.write(line)
            except IOError:
                print("error writing to revcomp fasta", IOError)
        else:
            for character in line:
                if character == "A":
                    line_string += "T"
                elif character == "T":
                    line_string += "A"
                elif character == "C":
                    line_string += "G"
                elif character == "G":
                    line_string += "C"
    try:
        with open("revcomp_fasta.fasta", "a") as revcomp_handle:
            if line_string != "":
                line_lengths.append(len(line_string)+prev_line_length)
                revcomp_handle.write(line_string[::-1])
    except IOError:
        print("error writing to revcomp_fasta", IOError)
    return line_lengths

# takes opened fasta file, number of non-stop codons before its an ORF, the species/ chromosome name specified by the user, the reading frame and the gff_lines list 
# detects open reading frames and calls write_to_fasta() function for each orf detected, writing both the header and the sequence of the ORF
# important information for the gff file for each ORF is stored in a list called gff_lines
def gene_search(fasta_handle, count, seqname, reading_frame, gff_lines, line_lengths=None):
    fasta_handle.seek(0)
    if(reading_frame > 2):
        reading_frame -= 3
        is_reverse = True
    else:
        is_reverse = False
    potential_gene = ""
    current_codon = ""
    genes_found = 0
    end_position = reading_frame
    waiting = 0
    header_num = -1
    for line in fasta_handle:
        # we can more or less ignore header lines
        if line[0] == ">":
            current_header = line
            waiting = 0
            header_num += 1
        else:
            for character in line:
                # check to see if we have waited enough characters for a given reading frame
                if(waiting < reading_frame):
                    waiting += 1
                elif(character != "\n"):
                    end_position += 1
                    current_codon += character
                    # if it is a 3 characters long, we check if this 'codon' is a stop codon
                    if len(current_codon) == 3:
                        potential_gene += current_codon
                        # if it is a stop codon and hits our confidence count, we record it as an ORF
                        if (current_codon == "TAA" or current_codon == "TAG" or current_codon == "TGA") and len(potential_gene) / 3 >= count:
                            start_position = end_position - len(potential_gene)
                            if(is_reverse):
                                gff_lines.append(GffLine(start=line_lengths[header_num]-start_position, end=line_lengths[header_num]-end_position, seqname=seqname, strand="-"))
                            else:
                                gff_lines.append(GffLine(start=start_position, end=end_position, seqname=seqname, strand="+"))
                            write_to_fasta(current_header, potential_gene)
                            potential_gene = ""
                            genes_found += 1
                        # if it is a stop codon that doesn't hit our confidence count we toss it
                        elif (current_codon == "TAA" or current_codon == "TAG" or current_codon == "TGA"):
                            potential_gene = ""
                        current_codon = ""

# writes potential genes (ORFs) to a file called blast_fasta.fasta which will be the input file for blast
def write_to_fasta(header, gene):
    try:
        with open("blast_fasta.fasta", "a") as blast_fasta:
            blast_fasta.write(header)
            blast_fasta.write(gene+"\n")
    except IOError:
        print("error writing to blast_fasta", IOError)

# runs the subprocess blast command on the command line with the file blast_fasta.fasta
def blast():
    print("starting blast")
    subprocess.run("blastn -query blast_fasta.fasta -db /home/share/databases/nt_db/nt -out blast_output_unfiltered.out -num_threads 10 -outfmt \"6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -max_hsps 1 -max_target_seqs 1", shell=True)
    #uncomment this and comment out line above to run blastx
    #subprocess.run("blastx -query blast_fasta.fasta -db /home/genome/shared/nr_db/nr -out blast_output_unfiltered.out -num_threads 10 -outfmt \"6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -max_hsps 1 -max_target_seqs 1", shell=True)
    print("finishd blast")

# parses the blast output file and pulls out the important information for the GFF file for each ORF and saves it into the list of GffLine objects
def parse_blast_output(gff_lines):
    try:
        with open("blast_output_unfiltered.out", "r") as blast_handle:
            key = 0
            for line in blast_handle:
                blast_line = line.split("\t")
                gff_lines[key].set_score(blast_line[3])
                gff_lines[key].set_attribute("{0}={1}".format(key, blast_line[2]))
                key += 1
    except IOError:
        print("error reading blast_output_unfiltered", IOError)

# removes any file using subprocess rm on the command line. mainly used for removing files at the beginning of the program that may have been created in previous runs
def remove_file(filename):
    subprocess.run("rm {0}".format(filename), shell=True)

# writes GffLine objects in gff_lines list to output file gff_output.gff in GFF format
def gff_writer(gff_lines):
    print("starting writing gff file")
    try:
        with open("gff_output.gff", "w") as gff_handle:
            for gff_line in gff_lines:
                gff_handle.write(str(gff_line))
                gff_handle.write("\n")
    except IOError:
        print("error writing to GFF", IOError)
    print("finished writing gff file")

# initializes class of GffLine which includes all 9 of the items that would be in a line of a gff file
class GffLine:
    def __init__(self, seqname=".", source="protein homology", feature="CDS", start=".", end=".", score=".", strand=".", phase="0", attribute="."):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attribute = attribute

    def set_seqname(self, seqname):
        self.seqname = seqname

    def set_source(self, source):
        self.source = source

    def set_feature(self, feature):
        self.feature = feature

    def set_start(self, start):
        self.start = start

    def set_end(self, end):
        self.end = end

    def set_score(self, score):
        self.score = score

    def set_strand(self, strand):
        self.strand = strand

    def set_phase(self, phase):
        self.phase = phase

    def set_attribute(self, attribute):
        self.attribute = attribute

    # returns a string version that fits gff file format
    def __str__(self):
        return self.seqname + "\t" + self.source + "\t" + self.feature + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + self.score + "\t" + self.strand + "\t" + self.phase + "\t" + self.attribute

# parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument("fasta_assembly", help="path to metagenomic assembly in fasta format", type=str)
parser.add_argument("seqname", help="name of the given sequence, usually either species or chromosome name", type=str)
parser.add_argument("--confidence_interval", help="desired confidence interval for gene detection", type=float, default=0.95) 

args = parser.parse_args()

# calling master function
annotate_genome(args.fasta_assembly, args.confidence_interval, args.seqname)

