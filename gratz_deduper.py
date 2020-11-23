#!/usr/bin/env python
#Derrik Gratz Fall 2020
import argparse
from os import remove
import subprocess
import re


def get_args():
    '''
    Store information provided by user at command line
    sample input: ./deduper.py -i $IN_DIR -o $OUT_DIR -u $UMI_FILE
    output:
        input_SAM = $IN_DIR
        output_directory = $OUT_DIR
        expected_umis = $UMI_FILE
    '''
    test_directory = '/home/dgratz/bgmp/deduper/Deduper/'
    parser = argparse.ArgumentParser(description='Removes PCR duplicates from a SAM file.\n\
        Duplicates are defined by having the same source group/chromosome, same starting position (adjusted for soft-clipping in the alignment),\
        same strand (forward or reverse), same unique molecular identifier (UMI), and same read length.\n\
        By default, the first encountered read is kept.\n This program does not currently support paired reads in the input SAM.\n\
        If you want to check for expected UMIs in reads, provide a list of UMIs in a single field text file with one UMI per line and nothing else.\n\
        Expected UMI checking is not necessary, and duplicate checking will still look for identical barcodes.')
    parser.add_argument('-u', '--umi', type=str, help='Absolute path of file containing expected UMIs', default=None)
    parser.add_argument('-i', '--input', type=str, help='Absolute path of SAM file of aligned reads', default=test_directory + 'test.sam')
    parser.add_argument('-o', '--output', type=str, help='Absolute path of directory to store output', default=test_directory + 'deduper_output')
    #add optional arg to change how duplicates are handled?
    parser.add_argument('-n', '--name',type=str, help='Optional prefix to name output sam files', default=None)
    parser.add_argument('-p', '--paired', help="Flag to denote input sam has paired reads", action="store_true")
    return parser.parse_args()


def sam_manipulation(input_sam, output_dir):
    '''
    Sorts the input sam file by read length, chromosome, then bitflag.
    Produces an intermediate sorted sam file that is used for duplicate checking, but is later deleted
    '''
    print("sorting the input file, this could take a while")
    sorted_file = output_dir + "sorted_input.sam"
    #sorting is performed with basic unix commands
    #adds the length of the read to the beginning of the file to facilitate sorting. This is later removed
    command = "cat {} | awk '{{print length($10)".format(input_sam) +'"\t"'+"$0;}' | unexpand -a" + '| sort -n -k 1,1 -k 4,4 -k 3,3 | cut -f 1 --complement | grep -v "^@" > {}'.format(sorted_file)
    subprocess.call(command, shell=True)
    print("Sorting complete!")
    return sorted_file


umi_checking = True
def grab_umis(umi):
    '''
    Makes a list of expected umis from the accessory file provided at runtime. If no file is provided, the program will not
    check for misindexed reads. It will still check UMIs when detecting duplicates. 
    Input: umi file
    output: a list of the expected UMIs
    '''
    global umi_checking
    umis = []
    if umi == None:
        print("No UMI file provided, skipping UMI checking")
        umi_checking = False
        return None
    with open(umi, "r") as fh:
        for line in fh:
            umis.append(line.strip())
    return umis


def cigar_adjustments(cigar, pos, strand):
    '''
    Calculates where the read would map to the reference if soft clipping consumed the reference
    Sample input: the starting position and cigar string from a single sam file read
        e.g. (1573, '15S75M', forward)
    Sample output: the adjusted starting position 
        e.g. 1558
    '''
    corrected_starting_pos = pos
    #split cigar string by number/character pieces
    stogies = re.findall('(\d+\w)', cigar)
    try:
        if strand == 'forward':
            if 'S' in stogies[0]:
                #Leading S indicates soft clipping, adjust by the number attached to S
                offset = int(re.search('\d+', stogies[0])[0])
                corrected_starting_pos = pos - offset
        elif strand == 'reverse':
            #reverse strand has inverted cigar string, need to look at end of cigar
            if 'S' in stogies[-1]:
                offset = int(re.search('\d+', stogies[-1])[0])
                corrected_starting_pos = pos - offset
        #might be able to remove this if can guarantee bitflags are set
        elif strand == 'undefined':
            #treating it as forward read for soft clipping purposes
            if 'S' in stogies[0]:
                offset = int(re.search('\d+', stogies[0])[0])
                corrected_starting_pos = pos - offset
    except IndexError:
        #Error: Improper cigar string. Can't check for soft clipping. Default to provided position
        return pos
    return corrected_starting_pos


unique_reads = set()
current_seqlen = 0
current_chrom = 0
def line_parser(line, expected_umis):
    '''
    Breaks apart a SAM file line and extracts relevant info for PCR duplicate detection
    Sample input: 1 non-header line from a sam file (e.g. NS500451:154:HWKTMBGXX:1:11101:4191:1194:TAGCAAGG	0	2	76718924	36	71M	*	0...)
    Sample output: whether it's a duplicated read, a unique read, or misindexed
    '''
    global current_seqlen
    global current_chrom
    global unique_reads
    global umi_checking
    try:
        fields = line.split('\t')
        qname = fields[0]
        umi = qname.split(':')[-1]
        #check if it's misindexed before doing more intense analyses
        if umi_checking:
            if umi not in expected_umis:
                return 'misindexed'
        bitflag = int(fields[1])
        #get strand from bitflag
        if ((bitflag & 16)==16):
            strand = 'forward'
        elif ((bitflag & 32)==32):
            strand = 'reverse'
        else:
            #may be able to remove 
            #print('Error: bitflag has neither forward nor reverse strand set')
            strand = 'undefined'
        chromosome = fields[2]
        if chromosome != current_chrom:
            #file was sorted by chromosome, a new chrom means new set of reads
            unique_reads.clear()
            current_chrom = chromosome
        pos = int(fields[3])
        cigar = fields[5]
        tlen = fields[8]
        seq_len = len(fields[9])
        if seq_len != current_seqlen:
            #the file has been sorted by the length of the read sequence, so if no sequences remain of that length, it can be assumed
            #there would be no more duplicates. The set of encountered reads can be reset. 
            unique_reads.clear()
            current_seqlen = seq_len
        adjusted_positions = cigar_adjustments(cigar, pos, strand)
        unique_identifiers = strand + str(chromosome) + str(adjusted_positions) + umi
        if unique_identifiers in unique_reads:
            return 'duplicate'
        else:
            unique_reads.add(unique_identifiers)
            return 'unique'
    except IndexError:
        #Something about this read is unexpected. Output to unqiue by default
        return 'unique'


def get_headers(input_sam):
    '''
    Get header lines from the input sam and write them to all the output files
    This saves a condition from being added to the loop that reads through the whole file
    '''
    headers = []
    checks = 0
    with open(input_sam, 'r') as fh:
        for line in fh:
            if line[0] == '@':
                headers.append(line)
            elif checks < 40:
                #check a few reads to make sure the file isn't paired end
                fields = line.split('\t')
                if ((int(fields[1]) % 1) == 1):
                    #first field of bitflag is true, read is paired
                    print("It looks like the input sam has paired reads. This program does not currently support paired-read SAM files. The program will now close.")
                    fh.close()
                    quit
                checks += 1
            else:
                fh.close()
                return headers


def main():
    '''
    Grab the UMIs from the file provided at runtime (if applicable)
    Open the input SAM passed in at runtime
    Sort the input SAM to produce an intermediate file
    Opens 3 output files, one for duplicate, misindexed, or unique
    Reads through the sorted input intermediate SAM file
    Passes the lines into the line parser to determine if the line is a PCR duplicate, misindexed, or unique
    write the line to the proper output file
    keep track of the lines added to each file
    '''
    global umi_checking
    args = get_args()
    #Not ready for paired reads yet
    if args.paired:
        print("It looks like the input sam has paired reads. This program does not currently support paired-read SAM files. The program will now close.")
        quit()
    #Grab expected umis from the command line (if applicable)
    expected_umis = grab_umis(args.umi)
    if expected_umis == None:
        umi_checking = False
    #making sure the directory ends in a slash
    if args.output[-1] != '/':
        args.output += '/'
    headers = get_headers(args.input)
    sorted_sam = sam_manipulation(args.input, args.output)
    #opening output files
    if args.name != None:
        output_prefix = args.output + args.name + "_"
    else:
        output_prefix = args.output
    outputs = {}
    outputs['misindexed'] = open((output_prefix + 'misindexed.sam'), "w")
    outputs['duplicate'] = open((output_prefix + 'duplicates.sam'), "w")
    outputs['unique'] = open((output_prefix + 'deduped.sam'), "w")
    #keeping track of how many reads are written to each file
    counts = {'misindexed':0, 'duplicate':0, 'unique':0}
    try:
        for line in headers: 
            for file in outputs:
                outputs[file].write(line)
    except TypeError:
        print("The header lines in the input file are not properly formatted. They will not be added to the output file.")
    #read through the sorted input file
    ln = 0
    with open(sorted_sam, 'r') as fh:
        print("Starting PCR duplicate checks")
        #no headers in the sorted file
        for line in fh:
            results = line_parser(line, expected_umis)
            outputs[results].write(line)
            counts[results] += 1
            ln += 1
            if ln % 1000000 == 0:
                print("Running: on read " + str(ln))
    for file in outputs:
        outputs[file].close()
    print('\nDuplication check complete! See summary below')
    for category in counts:
        print('Number of {} reads: {}   {:.2%}'.format(category, counts[category], counts[category]/ln))
    remove(sorted_sam)
    if umi_checking == False:
        remove(output_prefix + 'misindexed.sam')

main()
