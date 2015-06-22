#!/usr/bin/python

# import libraries
import sys # for argument vector
import optparse # to parse arguments

# describe what the script does
what_i_do = "A simple script to generate a ggbio genomics ranges file from an NCBI feature table."
usage = "myscript.py -i <input_file> -o <output_file>"

# initialize the parser
parser = optparse.OptionParser(usage = usage, epilog = what_i_do)
parser.add_option("-i", "--input_file", dest="input_file", default=None,
                   help='file to print out to the screen [Required]')
parser.add_option("-o", "--output_file", dest="output_file", default=None,
                  help='ggbio R genomic ranges data frame  [Required]')

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

# the main function of the script
def main():
    (opts, args) = parser.parse_args()
    # opts is a dictionary that contains all the options
    if not (opts.input_file and opts.output_file):
        print "Error: Input our output file not specified"
        print usage
        exit()

    # open the file 
    fh = open(opts.input_file, "r")
    fh_out = open(opts.output_file, "w")

    # read lines from file
    lines = fh.readlines()
    fh.close()
    
    header = ["chr", "start", "end", "strand", "itr_id", "gene", "product", "protein_id"]
    fh_out.write("\t".join(header) + "\n")
    
    # iterate through each line and print out content
    record = {} # record to fill in
    itr_id = 0 # unique record id
    flag = False
    for line in lines:
        line = line.strip("\n")
        fields = line.split("\t")
        if len(fields) == 3:
            if fields[2] == "CDS":
                if record and flag:
                    write_record(record, fh_out)
                    record = {}
                fr = fields[0]
                to = fields[1]
                cds = fields[2]
                strand = "+"
                ec="NA"
                if int(fr) > int(to):
                    strand = "-"
                    to = fields[0]
                    fr = fields[1]
                record["from"] = fr
                record["to"] = to
                record["cds"] = "CDS"
                record["strand"] = strand
                record["ec"] = ec
                record["itr_id"] = itr_id
                itr_id += 1
                if not flag:
                    # on first round skip the first gene
                    flag = True
        elif len(fields) > 3:
            if fields[3] == "gene":
                gene = fields[4]
                record["gene"] = gene
            if fields[3] == "product":
                product = fields[4]
                record["product"] = product
            if fields[3] == "EC_number":
                ec = fields[4]
                record["ec"] = ec
            if fields[3] == "protein_id":
                protein_id = fields[4]
                record["protein_id"] = protein_id
        else:
            continue
    
    
    fh_out.close()

if __name__ == "__main__":
    main()
