#!/usr/bin/env python
import argparse
#import pysam
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
    parser = argparse.ArgumentParser(description='Removes PCR duplicates from a SAM file of uniquely aligned reads. By default, the first encountered read is kept.')
    parser.add_argument('-u', '--umi', type=str, help='file containing expected UMIs', default='STL96.txt')
    parser.add_argument('-i', '--input', type=str, help='SAM file of aligned reads', default='test.sam')
    parser.add_argument('-o', '--output', type=str, help='Directory to store output', default='test_output')
    #add optional arg to change how duplicates are handled
    #add optional flags to correct reverse strand clipping 
    return parser.parse_args()


def sam_manipulation(input_sam, output_dir):
    sorted_file = output_dir + "sorted_input.sam"
    #subprocess.call('cat', input_sam, '| grep -v "^@" | sort -n -k 3,3 -k 9,9 -k 2,2 > sorted_input.sam')
    #
    #pysam.sort('-o', sorted_file, input_sam)
    
    return sorted_file


def grab_umis(umi):
    '''
    Makes a list of expected umis from the accessory file provided at runtime. If no file is provided, the program will not
    check for misindexed reads. It will still check UMIs when detecting duplicates. 
    Input: umi file
    output: a list of the expected UMIs
    '''
    umis = []
    with open(umi, "r") as fh:
        for line in fh:
            umis.append(line.strip())
    return umis


def cigar_adjustments(cigar, pos, strand):
    '''
    Calculates where the read would map to the reference if soft clipping and insertions consumed the reference
    Sample input: the starting position and cigar string from a single sam file read
        e.g. (1573, '15S75M', forward)
    Sample output: the adjusted starting position 
        e.g. 1558
    '''
    stogies = re.findall('(\d+\w)', cigar)
    if strand == 'forward':
        if 'S' in stogies[0]:
            offset = re.match('\d+', stogies[0])
            corrected_starting_pos = pos - offset
    elif strand == 'reverse':
        if 'S' in stogies[-1]:
            offset = re.match('\d+', stogies[-1])
            corrected_starting_pos = pos - offset
    #might be able to remove this if can guarantee bitflags are set
    if strand == 'undefined':
        corrected_starting_pos = pos
    return corrected_starting_pos


unique_reads = set()

current_tlen = '0'
stored_reads = {}
def line_parser(line):
    '''
    Breaks apart a SAM file line and extracts relevant info for PCR duplicate detection
    Sample input: 1 non-header line from a sam file (e.g. NS500451:154:HWKTMBGXX:1:11101:4191:1194:TAGCAAGG	0	2	76718924	36	71M	*	0...)
    Sample output: whether it's a duplicated read, a unique read, or misindexed
    '''
    global current_tlen
    fields = line.split('\t')
    print(fields)
    qname = fields[0]
    umi = qname.split(':')[-1]
    #check if it's misindexed before doing more intense analyses
    if umi not in expected_umis:
        return 'misindexed'
    bitflag = int(fields[1])
    #get strand from bitflag
    if ((bitflag & 16)==16):
        strand = 'Forward'
    elif ((bitflag & 32)==32):
        strand = 'Reverse'
    else:
        #may be able to remove 
        print('Error: bitflag has neither forward nor reverse strand set')
        strand = 'undefined'
    chromosome = fields[2]
    pos = int(fields[3])
    #pass in read.cigar instead? If line is from pysam
    #cigar = line.cigar
    cigar = fields[5]
    tlen = fields[8]
    seq_len = len(fields[9])
    if tlen != current_tlen:
        stored_reads.clear()
        current_tlen = tlen
    adjusted_positions = cigar_adjustments(cigar, pos, strand)
    unique_identifiers = 
    if unique_identifiers in unqiue_reads:
        return 'duplicate'
    else:
        unique_reads.add(unique_identifiers)
    return adjusted_positions


def main():
    args = get_args()
    #making sure the directory ends in a slash
    if args.output[-1] != '/':
        args.output += '/'
    sorted_sam = sam_manipulation(args.input, args.output)
    umis = grab_umis(args.umi)
    with open(args.input, 'r') as fh:
        for line in fh:
            if line[0] != '@':
                results = line_parser(line)
                print(line)

    #with pysam.AlignmentFile(sorted_sam, "r") as fh:
    #    file1 = 1

main()