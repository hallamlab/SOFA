#!/usr/bin/python

# import libraries
import sys # for argument vector
import optparse # to parse arguments
import re
import operator

# describe what the script does
what_i_do = "A simple script to generate a ggbio genomics ranges file from an NCBI feature table."
usage = "myscript.py -i <input_file> -o <output_file>"

# initialize the parser
parser = optparse.OptionParser(usage = usage, epilog = what_i_do)
parser.add_option("-r", "--read_sequences", dest="read_sequences", default=None,
                   help='File containing simulated reads from SOFAMetaSimCOG.py')
parser.add_option("-g", "--cds_genomic_ranges", dest="cds_genomic_ranges", default=None,
                   help='File containing genomic ranges of genome')
parser.add_option("-m", "--min_overlap", dest="overlap", default=0,
                  help='Amout of overlap between read and gene to be considered duplication [default=0]')
parser.add_option("-o", "--output_file", dest="output_file", default=None,
                  help='ggbio R genomic ranges with duplicated status')

def write_record(record, fh_out):
    line = ["NC_000913.3"]
    line.append(record["from"])
    line.append(record["to"])
    line.append(record["strand"])
    line.append(record["itr_id"])
    line.append(record["gene"])
    if "product" in record:
        line.append(record["product"])
    else:
        line.append("NA")
    if "protein_id" in record:
        line.append(record["protein_id"])
    else:
        line.append("NA")
    line = map(str, line) # change to string
    fh_out.write("\t".join(line) + "\n")


def process_cds_ranges(cds_genomic_ranges):
    cds_ranges = {}
    with open(cds_genomic_ranges) as gr_fh:
        header = gr_fh.readline()
        for line in gr_fh:
            line = line.strip("\n")
            fields = line.split("\t")
            cds_id = fields[4]
            fr = fields[1]
            to = fields[2]
            strand = fields[3]
            if cds_id not in cds_ranges:
                cds_ranges[cds_id] = {"from": fr, "to": to, "strand": strand}
    return cds_ranges

def process_read_sequences(read_sequences):
    read_pairs = {}
    with open(read_sequences) as rs_fh:
        pair_read_id = 0
        for line in rs_fh:
            if line.startswith(">"):
                results = SIM_PATTERN.search(line)
                if results:
                    pair_id = results.group(1)
                    fr = results.group(2)
                    to = results.group(3)
                    pair_num = results.group(4)
                    
                    if pair_id not in read_pairs:
                        read_pairs[pair_id] = {}
                    if pair_num not in read_pairs[pair_id]:
                        read_pairs[pair_id][pair_num] = {}
                    read_pairs[pair_id][pair_num]["from"] = fr
                    read_pairs[pair_id][pair_num]["to"] = to
    return read_pairs

def is_overlap(read, cds, overlap=0):
    # Checks to see if there is an overlap between the read and the CDS
    cds_from = int(cds["from"])
    cds_to = int(cds["to"])
    read_from = int(read['from'])
    read_to = int(read['to'])
    
    if read_from < cds_from and read_to < cds_from:
        return False
    elif read_from > cds_to and read_to > cds_to:
        return False
    else:
        # print cds_from, cds_to, read_from, read_to, overlap
        if read_to > (cds_from + overlap) and read_from < cds_from:
            return True
        elif read_from < (cds_to - overlap) and read_to > read_from:
            return True
        elif read_from > cds_from and read_to < cds_to:
            return True
        else:
            return False
    

def is_double_counting(rp, cds, overlap):
    # checks to see if there is double counting in a paired-end read and
    # a particular CDS
    count = 0
    for read in rp:
        hitset = []
        if is_overlap(rp[read], cds, overlap=overlap):
            count += 1
    if count == 2:
        return True
    return False

def write_rp(rp_id, rp, cds_hits, out_fh):
    # writes out duplication genomic ranges
    chrome = "NC_000913.3"
    strand = "+" # assumed forward strand
    fr = rp['1']['from']
    to = rp['2']['to']
    cds_ids = "NA"
    dbl_count = "FALSE"
    if len(cds_hits) > 0:
        # take from of read one to the to of read two as genomic range coords
        cds_ids = ",".join(cds_hits)
        dbl_count = "TRUE"

    # write line
    line = [chrome, fr, to, strand, dbl_count, cds_ids]
    line = "\t".join(line) + "\n"
    out_fh.write(line)

# pattern for parsing reads
SIM_PATTERN = re.compile("##(.+?)_(.+?):(.+?)_([0-9]+?)$")

def genome_coverage(cds_ranges):
    # convert CDS ranges
    L = 4641652 # E.coli genome length
    
    my_items = cds_ranges.items()
    temp_key_start_coords = list()
    for item in my_items:
        temp_key_start_coords.append((int(item[0]), int(item[1]['from'])))
    
    # sort by starting position
    cds_ranges_sorted = sorted(temp_key_start_coords, key=operator.itemgetter(1))
    
    c = 0
    l = 0 # genome length
    for itr in cds_ranges_sorted:
        a = int(cds_ranges[str(itr[0])]["from"])
        b = int(cds_ranges[str(itr[0])]["to"])
        if (a-c) > 0:
            l += (a-c)
        c = b
    l += (L-c)
    print "Non-coding DNA:", l
    print "Coding DNA:", L-l
    print "Coding Density:", float(L-l) / float(L)
    

# the main function of the script
def main():
    (opts, args) = parser.parse_args()
    # opts is a dictionary that contains all the options
    if not (opts.read_sequences and opts.cds_genomic_ranges and opts.output_file):
        print "Error: Input our output file not specified"
        print usage
        exit()
    
    # parse read sequences
    read_pairs = process_read_sequences(opts.read_sequences)

    # parse CDS ranges
    cds_ranges = process_cds_ranges(opts.cds_genomic_ranges)
    genome_coverage(cds_ranges)
    
    # header
    header = ["chr", "start", "end", "strand", "dbl_count", "cds"]
    
    out_fh = open(opts.output_file, 'w')
    out_fh.write("\t".join(header) + "\n")
    
    # figure out which read pairs are double counting
    for rp in read_pairs:
        cds_hits = []
        for cds in cds_ranges:
            if is_double_counting(read_pairs[rp], cds_ranges[cds], int(opts.overlap)):
                cds_hits.append(cds) # add cds_ids to list
                
        write_rp(rp, read_pairs[rp], cds_hits, out_fh)
    
    
    exit()

if __name__ == "__main__":
    main()
