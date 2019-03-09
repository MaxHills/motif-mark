#!/usr/bin/env python3
# Python3 must be loaded prior to implementing this program.
#######################
## IMPORT STATEMENTS ##
#######################

import argparse
from argparse import RawTextHelpFormatter
import re
import cairocffi as cairo
import math
import random

##########################
## FUNCTION DEFINITIONS ##
##########################

def readMotifs(file):
    """
    Input: A text file with 1 motif sequence per line
    Output: Returns a tuple with: an int number of motifs, a list of motif sequence strings,
    a list of integers representing the nucleotide length of each motif, and an re search
    term, per motif
    """
    motifs = []
    with open(file, 'r') as f:
        re_terms = []
        for line in f:
            line = line.strip()
            line = line.upper()
            this_term = ""
            for char in line:
                # Create an re search term, this_term, for each motif
                if char in "ACTG":
                    this_term += char
                elif char == "U":
                    this_term += "T"
                elif char == "Y":
                    this_term += "[CT]"
            motifs.append(line)
            re_terms.append(re.compile(this_term))
    num_motif = len(motifs)
    return num_motif,motifs,re_terms

def readSeqs(file):
    # Initialize dictionary for output
    seq_dict = {}
    # Initialize a list for sequences
    sequences = []
    # A 0 at pos 0 of list, seq_start, indicates sequence 0 starts with an intron.
    # A 1 at pos 0 of seq_start indicates sequence 0 starts with an exon.
    seq_start = []
    # The lengths of each sequence
    seq_lens = []
    # A list of intron lengths for each sequence
    intr_lens = []
    # A list of exon lengths for each sequence
    exon_lens = []
    # A list of intron start positions for each seq
    intr_strt = []
    # A list of exon start positions for each seq
    exon_strt = []
    # A list of gene names
    gene_names = []
    with open(file, 'r') as f:
        LN = 0
        for line in f:
            # Special case for first line of file
            if LN == 0:
                LN += 1
                # Grab the gene name
                gene_names.append(re.findall(">([A-Z0-9]+)\s", line))
                # Initialize Sequence
                this_seq = ""
            # The next gene & sequence
            elif line.startswith(">"):
                # Add this_seq to the list, sequences
                sequences.append(this_seq)
                # Grab the gene name
                gene_names.append(re.findall(">([A-Z0-9]+)\s", line))
                # Initialize Sequence
                this_seq = ""
            else:
                # Remove newline characters
                line = line.strip()
                # Add the line to the current sequence
                this_seq = this_seq + line
        # Add the last sequence in the file to sequences
        sequences.append(this_seq)
        # How many sequences?
        num_seqs = len(gene_names)
        for i in range(num_seqs):
            this_in_lens = []
            this_ex_lens = []
            this_in_pos = []
            this_ex_pos = []
            # Does the sequence begin with an INTRON (0) or an EXON (1)?
            if sequences[i][0] in "actg":
                seq_start.append(0)
            else:
                seq_start.append(1)
            seq_lens.append(len(sequences[i]))
            intr_srch = re.compile("[actg]+")
            introns = intr_srch.finditer(sequences[i])
            for intron in introns:
                this_in_lens.append(len(intron.group()))
                this_in_pos.append(intron.start())
            exon_srch = re.compile("[ACTG]+")
            exons = exon_srch.finditer(sequences[i])
            for exon in exons:
                this_ex_lens.append(len(exon.group()))
                this_ex_pos.append(exon.start())
            intr_lens.append(this_in_lens)
            intr_strt.append(this_in_pos)
            exon_lens.append(this_ex_lens)
            exon_strt.append(this_ex_pos)
        print("Raw INTRON STARTS:\n",intr_strt)
    for n in range(num_seqs):
        this_gene = gene_names[n][0]
        seq_dict[this_gene] = (seq_lens[n],sequences[n],intr_strt[n],intr_lens[n],exon_strt[n],exon_lens[n],seq_start[n])
    return seq_dict

def findMotifs(gene_names, motif_seq, motif_reLIST, sequenceDICT):
    # A dictionary of dictionaries: KEY=Gene Name, VALUE=seq_pos Dictionaries, see seq_pos comment
    seq_mot_pos = {}
    for key in sequenceDICT:
        this_seq = sequenceDICT[key]
        # A dictionary of lists: KEY=Motif Sequence, VALUE=A list of Motif Starting Positions for each Sequence
        seq_pos = {}
        for motif in motif_reLIST:
            j = motif_reLIST.index(motif)
            # A list of starting positions for each motif
            mot_pos = []
            # An iterator object of integer starting positions
            positions = motif.finditer(this_seq.upper())
            for p in positions:
                mot_pos.append(p.start())
            motif_key = motif_seq[j]
            seq_pos[motif_key] = mot_pos
        seq_mot_pos[key] = seq_pos
    return seq_mot_pos
    
def rgbGEN(n):
    random.seed(13)
    rgb_list = []
    for _ in range(n):
        rgb = [random.uniform(0,1) for _ in range(3)]
        rgb_list.append(rgb)
    return rgb_list

def rgbGENex(n):
    random.seed(12354734)
    rgb_list = []
    for _ in range(n):
        rgb = [random.uniform(0,1) for _ in range(3)]
        rgb_list.append(rgb)
    return rgb_list
    
def drawIntron(context,x_1,y,x_2):
    context.set_line_width(10)
    context.set_source_rgb(0,0,0)
    context.move_to(x_1,y)
    context.line_to(x_2,y)
    context.stroke()
    
def drawExon(context,exon_col,x_1,y,x_2):
    # Motif-Top
    y_1 = y - 10
    # Motif-Bottom
    y_2 = y + 10
    context.set_source_rgb(exon_col[0],exon_col[1],exon_col[2])
    context.set_line_width(2)
    print("Exon Coordinates:\n"+str(x_1)+"\n"+str(y_1)+"\n"+str(x_2)+"\n"+str(y_2))
    context.move_to(x_1,y_1)
    context.line_to(x_2,y_1)
    context.line_to(x_2,y_2)
    context.line_to(x_1,y_2)
    context.line_to(x_1,y_1)
    context.fill()
    context.stroke()
    
def drawMotif(context,motif,positionLIST,motif_color,y):
    bottom = y + 10
    top = y - 10
    context.set_source_rgb(motif_color[0],motif_color[1],motif_color[2])
    # Set the line-width to the width of the motif.
    linwid = len(motif)
    # Shift the motif half its width, to center.
    shift = linwid//2
    context.set_line_width(linwid)
    for position in positionLIST:
        context.move_to(position+shift,top)
        context.line_to(position+shift,bottom)
        context.stroke()
    
def drawGene(context,gene_name,intr_strts,intr_lens,exon_strts,exon_lens,exon_color,
            mot_posDICT,motif_colors,y, width, gene_length):
    # posit is a dictionary: KEY=MotifSeq, VALUE=MotifPositions
    x_shift = (width-gene_length)//2
    motif_names = list(mot_posDICT.keys())
    for name in motif_names:
        if len(mot_posDICT[name]) > 0:
            n = len(mot_posDICT[name])
            # pos is a list of positions
            for pos in range(n):
                mot_posDICT[name][pos] += x_shift
########################################################
########################################################
    # Create Intron X_2 values
    intr_x2s = []
    for intr in range(len(intr_strts)):
        intr_x2s.append(intr_lens[intr] + intr_strts[intr])
    # Create Exon X_2 values
    exon_x2s = []
    for exon in range(len(exon_strts)):
        exon_x2s.append(exon_lens[exon] + exon_strts[exon])
    # Draw the introns
    context.translate(0,100)
    for intr in range(len(intr_strts)):
        # drawIntron(context,x_1,y,x_2)
        print("Intron X1s:\n"+str(intr_strts)+"\nIntrons X2s:\n"+str(intr_x2s))
        drawIntron(context,intr_strts[intr],y,intr_x2s[intr])
    # Draw the exons
    exon_col = exon_color[gene_name]
    for exon in range(len(exon_strts)):
        # drawExon(context,exon_color,x_1,y,x_2)
        drawExon(context,exon_col,exon_strts[exon],y,exon_x2s[exon])
    # Draw all motifs, per gene
    for mot in motif_names:
        pos_list = mot_posDICT[mot]
        mot_col = motif_colors[mot]
        # drawMotif(context,motif,positionLIST,motif_color,y)
        drawMotif(context,mot,pos_list,mot_col,y)
    # Determine gene-sequence center
    xtxtcntr = x_shift + (gene_length//2)
    # Write the gene name
    context.select_font_face('Sans')
    context.set_font_size(15) # em-square height is 90 pixels
    context.move_to(xtxtcntr-15,y-20) 
    context.set_source_rgb(0, 0, 0) # black
    context.show_text(gene_name)

#####################
## ARGUMENT PARSER ##
#####################

# Create the argument parser.
parser = argparse.ArgumentParser(
        description=
        """
        INPUT:  A text file with one motif per line
                A FASTA file of sequences to search for the given motifs;
                Intronic sequence should be lowercase and exonic in UPPERCASE
        OUTPUT: A scaled PDF or SVG image showing sequences with introns,
                exons, and motif locations
        """, formatter_class = RawTextHelpFormatter)
# Inform the parser of command line options.
# -h,--help is an argparse built-in function and will return the description w/ 
# options and option help messages seen below.
parser.add_argument('-f', '--fasta', help='Absolute FASTA file path', 
                    type=str, required=True)
parser.add_argument('-m', '--motif', help='Absolute motif file path; with one motif per line', 
                    type=str, required=True)
parser.add_argument('-p', '--pdf', help='Designates output image should be in PDF format',
                    action='store_true', required=False)
# Call the parser.
args = parser.parse_args()
# Make sure files were provided by user, if not, exit program with error message.
if args.fasta == None:
    msg = """\nERROR: A FASTA file of sequences must be provided
    To specify your FASTA filename use: -f,--fasta <FASTAfilename>
    \n"""
    print(msg)
    exit()
if args.motif == None:
    msg = """\nERROR: A motif sequence OR a file of motif sequences must be provided
    To specify your motif filename use: -m,--motif <MOTIFfilename>
    \n"""
    print(msg)
    exit()
# A tuple
motif_info = readMotifs(args.motif)
# A dictionary
seq_info = readSeqs(args.fasta)
# Extract sequences from seq_info
# and add gene names to a list
sequences = {}
geneNames = []
for key in seq_info:
    sequences[key] = seq_info[key][1]
    geneNames.append(key)
# geneNames = list of gene names, motif_info[1] = list of motif sequences, motif_info[3] = list of motif search terms
# sequences = A dictionary, keys=gene names, values=gene sequences
positions = findMotifs(geneNames,motif_info[1],motif_info[2],sequences)

### VARIABLES ##
motif_num = motif_info[0]
motif_names = motif_info[1]
motif_search = motif_info[2]

gene_names = list(seq_info.keys())
for name in gene_names:
    
    i = gene_names.index(name)
    

##############################################
##  ->->LET'S DRAW   #########################
##############################################

##  DETERMINE WIDTH & HEIGHT  ##
# Hardcode Width based on largest sequence length and screen limits.
width = 800
# Determine appropriate canvas height
seq_num = len(seq_info)
height = ((motif_num*50) + (seq_num*100))

##  SET THE SURFACE & CONTEXT  ##
if args.pdf == True:
    surface = cairo.PDFSurface("motif-mark.pdf",width,height)
    context = cairo.Context(surface)
else:
    surface = cairo.SVGSurface("motif-mark.svg",width,height)
    context = cairo.Context(surface)

## Encircle the canvas, to see dimensions
context.set_line_width(5)
context.set_source_rgb(0,0,0)
context.rectangle(0,0,width,height)
context.stroke()
    
##  DRAW A LEGEND  ##
x,y = 5,5
motif_colors = {}
for name in motif_names:
    i = motif_names.index(name)
    motif_colors[name] = rgbGEN(motif_num)[i]
# Draw Colored Squares
for key in motif_colors:
    context.set_line_width(1)
    the_cols = motif_colors[key]
    context.set_source_rgb(the_cols[0],the_cols[1],the_cols[2])
    context.rectangle(x,y,(x+15),(y+15))
    context.fill()
    # Advance x for Motif Sequence Text
    # draw text
    context.select_font_face('Sans')
    context.set_font_size(15) # em-square height is 90 pixels
    context.move_to(x+25,y+15) # move to point (x, y) = (10, 90)
    context.set_source_rgb(0, 0, 0) # black
    context.show_text(key)
    # Shift origin for next gene
    context.translate(0,25)

y = context.get_current_point()[1]

# Populate Exon Colors    
exon_colors = {}
exon_names = list(seq_info.keys())
for gene in exon_names:
    i = exon_names.index(gene)
    exon_colors[gene] = rgbGENex(seq_num)[i]

# PARAMETERS for drawGene:
# gene_name,gene_len,intr_strt,intron_lenLIST,exon_strt,exon_lenLIST,
# start_type,window_width,y,context,pos_dict
for key in seq_info:
    curx, cury = context.get_current_point()
    print("##########\nCurrent X,Y:\n"+str(curx)+","+str(cury))
    geneLEN = seq_info[key][0]
    intrSTRT = seq_info[key][2]
    print("In Main Program INT STRT:\n",str(intrSTRT))
    intronLENLIST = seq_info[key][3]
    exonSTRT = seq_info[key][4]
    exonLENLIST = seq_info[key][5]
    STRT_type = seq_info[key][6]
    x_shift = (width - geneLEN)//2
####### Shift all start positions based on x_shift #####
########################################################
    for start in range(len(intrSTRT)):
        intrSTRT[start] += x_shift
    for start in range(len(exonSTRT)):
        exonSTRT[start] += x_shift
    # pos_shifted is a dictionary: KEY=MotifSeq, VALUE=MotifPositions
    pos_shifted = positions[key]
    # drawGene(context,gene_name,intr_strts,intr_lens,exon_strts,exon_lens,exon_color,
             # positionDICT,motif_colors,y, width, gene_length)
    drawGene(context,key,intrSTRT,intronLENLIST,exonSTRT,exonLENLIST,exon_colors,
             pos_shifted,motif_colors,y,width,geneLEN)
    
    
surface.finish()    
    
    


print("Motif Information:\n",motif_info)
print("Sequence Information:\n",seq_info)
print("Positions:\n",positions)
print("Shifted Positions:\n", str(pos_shifted))
print("Motif Colors:", motif_colors)


