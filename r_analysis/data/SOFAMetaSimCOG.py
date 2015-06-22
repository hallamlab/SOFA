# Aria Hahn March 27 2015

import sys
import argparse
import re
from random import randint, seed

class FastaRecord():
    def __init__(self, longname, sequence):
      self.longname = longname
      self.sequence = sequence
      fields = [ x.strip() for x in self.longname.split(' ') ]
      if len(fields) > 0:
         self.name = fields[0]
      else:
         self.name = None

class FastaReader():
    """Parses a fasta record from a string or file."""
    stop = False
    START_PATTERN = re.compile(r'^>')
    name = None
    future_name =None
    sequence=""
    def __init__(self, fasta_filename):
        try:
            self.file = open(fasta_filename, 'r')
        except IOError:
            print "Cannot open fasta file " + fasta_filename

    def __iter__(self):
        return self

    def next(self):
        if self.stop:
          raise StopIteration
        
        try:
           if not self.name: 
               self.name = self.file.readline().strip()
           line = self.file.readline()
        except:
           line = None

        if not line:
           self.stop = True
           raise StopIteration

        fragments = []
        while line and not self.START_PATTERN.search(line):
            fragments.append(line.strip()) 
            line = self.file.readline()

       # print line
        if self.future_name:
            self.name = self.future_name

        if line:
            self.future_name = line.strip()
        
        self.sequence =''.join(fragments)
        self.seqname = self.name
        
        return FastaRecord(self.name, self.sequence)

def get_parser():
    parser = argparse.ArgumentParser(
        description='simulates sequencing from a fasta file given your parameters - Note that there is no error model and sequences shorter then the length of an input sequence - (insert length + range)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input_sequence', help='The name of the input fasta file')
    parser.add_argument('-o', '--output', help='The name of the output')
    parser.add_argument('-n', '--pairnum', help='the number of read pairs from EACH input sequnce in the fasta file')
    parser.add_argument('-l', '--readlen', help='the length of the reads') 
    parser.add_argument('-t', '--insertlen', help='the length of the pair reads - this can help you to determine if the reads overlap or not')
    parser.add_argument('-r', '--range', help='reads will be readlenth + (0-range). Set to zero if youd like uniform read lengths')
    return parser
    

def main():
    args = get_parser().parse_args()
    #print args.input_sequence
    
    fasta_reader = FastaReader(args.input_sequence)
    
    if args.output:
        fp = open(args.output, 'w')
    
    count = int(args.pairnum)
    readl = int(args.readlen)
    insert = int(args.insertlen)
    rangeby = int(args.range) 
    list = {}
    #print count
    #print readl
    # seed(12345) # first estimate
    # seed(54321) # second estimate
    seed(13232) # third estimate
    
    pair_counter = 0
    
    for record in fasta_reader:
        name = record.longname[1:]
        line = record.sequence
        count = int(args.pairnum)
        length=(len(line))#print length
        
        while count > 0:
            if (length-(insert+rangeby)) > 0:
                rand = randint(0,(length-(insert+rangeby)))
                rand2 = randint(0,rangeby)
                #print rand2    
                rand3 = randint(0,rangeby)
                #print rand3
                # print str(rand + readl + rand2)
                if rand not in list:
                    fp.write(">" + name + "##" + str(pair_counter) + "_" + str(rand) + ":" + str(rand + readl + rand2) + "_1" + "\n")
                    fp.write(line[(rand):(rand + readl + rand2)]+ "\n")
                    fp.write(">" + name + "##" + str(pair_counter) + "_" + str(rand + (insert -readl)) + ":" + str(rand + (insert -readl) + readl + rand3 ) + "_2" + "\n")
                    fp.write(line[(rand + (insert-readl)):(rand + (insert -readl) + readl + rand3)] + "\n")
                    pair_counter += 1
                    list[rand] = 1
                    count = count - 1
                # list = {}
            else:
                count = count -1
        list = {}
    # close the files
    fp.close()



if __name__ == '__main__':
    main()
