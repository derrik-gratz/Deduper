#!/usr/bin/env python
#Derrik Gratz Fall 2020
import argparse
import os
import pysam
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
        Expected UMI checking is not necessary, and duplicate checking will still look for identical barcodes.\n\
        A pre-sorted SAM can be passed in, or if the "-s" flag is not set, the program will attempt to sort the input SAM')
    parser.add_argument('-f', '--file', type=str, help='Absolute path of SAM file of aligned reads', required=True)
    parser.add_argument('-o', '--output', type=str, help='Absolute path of directory to store output', default=os.getcwd())
    parser.add_argument('-u', '--umi', type=str, help='Absolute path of file containing expected UMIs, optional for duplicate checking, but needed to check for misindexing', default=None)
    #add optional arg to change how duplicates are handled?
    parser.add_argument('-n', '--name', type=str, help='Optional prefix to name output sam files', default=None)
    parser.add_argument('-p', '--paired', help="Flag to denote input sam has paired reads, currently causes program to exit", action="store_true")
    parser.add_argument('-s', '--sorted', help="Flag to denote input sam is already sorted by chromosome", action="store_true")
    return parser.parse_args()


def sam_manipulation(input_sam, output_dir, sorting):
    '''
    Sorts the input sam file by read length, chromosome, then bitflag.
    Produces an intermediate sorted sam file that is used for duplicate checking, but is later deleted
    '''
    #check runtime flag first
    if not sorting:
        print("sorting the input file, this could take a while")
        sorted_file = output_dir + "sorted_input.sam"
        try:
            pysam.sort("-o", sorted_file, input_sam)
        except:
            #incase pysam isn't available on that computer
            print('An error occured with pysam sort. Resorting to UNIX sort')
            command = "cat {}".format(input_sam) + ' | unexpand -a | sort -n -k 3 | grep -v "^@" > {}'.format(sorted_file)
            subprocess.call(command, shell=True)
        print("Sorting complete!")
        return sorted_file
    else:
        print("Input file has been designated as sorted.")
        return input_sam


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
                corrected_starting_pos -= offset
        elif strand == 'reverse':
            #need to calculate the 5' starting position from the leftmost provided position
            #reverse strand has inverted cigar string, need to look at end of cigar for soft clipping
            if 'S' in stogies[-1]:
                offset = int(re.search('\d+', stogies[-1])[0])
                corrected_starting_pos += offset
            #looking for other consuming bases
            for stogie in stogies:
                if (re.search('\D', stogie)[0]) in 'DNM=X':
                #expanded funtionality for newer SAM cigars
                    offset = int(re.search('\d+', stogie)[0])
                    corrected_starting_pos += offset
                #currently unclear how best to implement this or if it's needed. Suspended functionality for now but leaving for posterity.
                #elif 'I' in stogie:
                #    offset = int(re.search('\d+', stogie)[0])
                #    corrected_starting_pos -= offset
    except IndexError:
        #Can't check for soft clipping. Default to provided position
        print('Error: unexpected cigar string {}. Cannot adjust starting position for soft clipping.'.format(cigar))
        return pos
    #if strand is undefined (neither forward nor reverse flag set), it just return starting position too
    return corrected_starting_pos


unique_reads = set()
current_chrom = str(0)
def line_parser(line, expected_umis):
    '''
    Breaks apart a SAM file line and extracts relevant info for PCR duplicate detection
    Sample input: 1 non-header line from a sam file (e.g. NS500451:154:HWKTMBGXX:1:11101:4191:1194:TAGCAAGG	0	2	76718924	36	71M	*	0...)
    Sample output: whether it's a duplicated read, a unique read, or misindexed
    '''
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
        if ((bitflag & 32)==32):
            strand = 'forward'
        elif ((bitflag & 16)==16):
            strand = 'reverse'
        else:
            #print('Error: bitflag has neither forward nor reverse strand set')
            strand = 'undefined'
        chromosome = fields[2]
        if chromosome != current_chrom:
            #file was sorted by chromosome, a new chrom means new set of reads
            unique_reads.clear()
            current_chrom = chromosome
        pos = int(fields[3])
        cigar = fields[5]
        adjusted_positions = cigar_adjustments(cigar, pos, strand)
        unique_identifiers = strand + str(chromosome) + str(adjusted_positions) + umi
        if unique_identifiers in unique_reads:
            return 'duplicate'
        else:
            unique_reads.add(unique_identifiers)
            return 'unique'
    except IndexError:
        print('Error in line, cannot check for uniqueness')
        print(line)
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
                if ((int(fields[1]) & 1) == 1):
                    #first field of bitflag is true, read is paired
                    print("\nIt looks like the input sam has paired reads. This program does not currently support paired-read SAM files. The program may not remove duplicates as would be expected with paried end reads, but it will attempt to remove duplicates while treating reads as single end.")
                    #quit()
                    break
                checks += 1
            else:
                break
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
        print("It looks like the paired flag is set. This program does not currently support paired-read SAM files.\n\
            If you wish to attempt deduplication on the file, you can pass it in again without the --paired flag.\n\
            This will treat the file as single end for the purposes of duplication checks.\n\
            This is not recommended, as the deduplication will not be accurate, but the program will still be functional.\n\
            The program will now close.")
        quit()
    #Grab expected umis from the command line (if applicable)
    expected_umis = grab_umis(args.umi)
    if expected_umis == None:
        umi_checking = False
    #making sure the directory ends in a slash
    if args.output[-1] != '/':
        args.output += '/'
    headers = get_headers(args.file)
    sorted_sam = sam_manipulation(args.file, args.output, args.sorted)
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
            if not line.startswith('@'):
                #results are whether the line is a duplicate, misindexed, or unique
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
        print('Number of {} reads:\t{}\t{:.2%}'.format(category, counts[category], counts[category]/ln))
    if not args.sorted:
        os.remove(sorted_sam)
    if umi_checking == False:
        os.remove(output_prefix + 'misindexed.sam')

main()
