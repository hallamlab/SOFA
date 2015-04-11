#!/usr/bin/python


try:
   import optparse, sys, re, csv, traceback, pickle
   from os import path, _exit, remove, rename, makedirs, popen 
   import logging.handlers
   from glob import glob
#   from libs.python_modules.utils.sysutil import pathDelim
   #from libs.python_modules.utils.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     print traceback.print_exc(10)
     sys.exit(3)

"""This script runs the SOFA pipeline """
parser = None
PATHDELIM = "/"

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


def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell."""
    pipe = popen(cmd + ' 2>&1', 'r')
    text = pipe.read()
    sts = pipe.close()
    if sts is None: sts = 0 
    if text[-1:] == '\n': text = text[:-1]
    return sts, text

def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def eprintf(fmt, *args):
    sys.stderr.write(fmt % args)

def files_exist( files , errorlogger = None):
    status = True    
    for file in files:
       if not path.exists(file):
          if errorlogger:
             errorlogger.write( 'ERROR\tCould not find ptools input  file : ' +  file )
          status = False
    return not status


def createParser():
    global parser

    usage = sys.argv[0] + """ -i <inputfolder> -s <sampleName> -o <outputFolder> --FlashExec <FLASH_executable>\n""" +\
                      """     --FragGeneScanExec <FragGeneScan_executable> --LASTExec <LAST_executable> --stage s \n(optional flags: --tempdirs t --bitScore bs  --evalue ev --tempdirs <tempdirs>\n"""

    epilog = """This script takes an interleved FASTQ file and uses the SOFA pipeline:

             Stage 1. : SOFA uses FLASH to merge pairs of reads when possible and produces two
			 FASTQ file containing merged and unmerged reads

             Stage 2. : SOFA concatenates the FASTQ files from stage 1 and converts
			 the FASTQ to FASTA format and provides a mapping file

             Stage 3. : Predicts ORFs using FragGeneScan+ and produces a .faa file 

             Stage 4. : The resulting unmerged reads in the .faa file are LASTED against a 
                        clustered reference protein database (REFSEQ proteins at 85% similarity). 
                        If both the reads of a pairs have hits with the same protein function then
                        one of the reads is removed. The final ORFs are in a .FINAL.faa file
             """

    #epilog = re.sub(r'\s+', ' ', epilog)

    #parser = optparse.OptionParser(usage=usage, epilog = epilog)
    parser = optparse.OptionParser(usage= usage + epilog)

    parser.add_option('-i', '--input', dest='inputfolder', default=None,
                           help='the path to the input folder')

    parser.add_option('-s', '--sample', dest='sample_name', default=None,
                           help='sample name - NO file extension')

    parser.add_option('-o', '--output', dest='outputfolder', default=None,
                           help='the output foldername')

    parser.add_option('--FlashExec', dest='FlashExecutable', default=None,
                           help='FLASH Executable')

    parser.add_option('--LASTExec', dest='LASTExecutable', default=None,
                           help='LAST Executable')

    parser.add_option('--refdb', dest='refdb', default=None,
                           help='RefDBs for the de-deduplication')

    parser.add_option('--FragGeneScanExec', dest='FragGeneScanExecutable', default=None,
                           help='FragGeneScan Executable')

    #parser.add_option('-M', dest='map_file', default=None, help='The map file')

    #parser.add_option('-p', dest='preprocessed_file', default=None, help='The preprocessed file')

    parser.add_option('--stage', dest='stage', action='append', default=[], choices=['1', '2', '3', '4' ], 
                           help='The stages to execute 1 : FLASH; 2 : Format files; 3 : FragGeneScan+ and  4 : Deduplication ')

    parser.add_option('--threads', dest='threads', default=1,
                           help='the number of threads to use with FragGeneScan+')
    
    parser.add_option('--bitScore', dest='bitScore', default=20,
                           help='Optional bit-Score cutoff value (default is 20)')
                           
    parser.add_option('--evalue', dest='evalue', default=0.000001,
                           help='Optional e-value cutoff value (0.000901)')
                        
    parser.add_option('--tempdirs', dest='tempdirs', action='append', default=[],  
                           help='the temp dirs to use - default is usually fine (no need to specify) ')



inputfiles = []

def  _execute_LAST(lastargs):
    args= [ ] 
    
    if lastargs['last_executable'] :
       args.append( lastargs['last_executable'])
	
    args += [ "-f",'2' ]
    
    if lastargs['bitScore']:
       args += [ "-S", lastargs['bitScore'] ]
       
    if lastargs['evalue']:
       args += [ "-E", lastargs['evalue'] ]
    
    if lastargs['last_output']:
       args += [ "-o", lastargs['last_output'] + ".tmp"]

    if lastargs['refdb']:
       args += [ lastargs['refdb'] ]
      
    if lastargs['query_file']:
       args += [ lastargs['query_file'] ]

    try:
       result = getstatusoutput(' '.join(args) )
       
       rename( lastargs['last_output']+ ".tmp",lastargs['last_output']) 
    except:
       message = "Could not run LASTAL correctly"
       if result and len(result) > 1:
          message = result[1]
       return (1, message)

    return (result[0], result[1])
      


def checkInput(options):
    global inputfiles
    #print options.inputfolder, options.sample_name, options.outputfolder, options.FlashExecutable

    errors = [] 
    if not path.exists(options.inputfolder):
       errors.append( "ERROR : Input folder " + options.inputfolder + " does not exist ")

    inputfile =  options.inputfolder + PATHDELIM + options.sample_name + ".fastq"

    if not inputfile:
       errors.append( "ERROR : Sample in folder  " + inputfile +  "  do not exist ")

    setattr(options, 'inputfile', inputfile)

    if not path.exists(options.FlashExecutable):
       errors.append( "ERROR : FLASH executable " + options.FlashExecutable + " does not exist ")

    if not path.exists(options.LASTExecutable):
       errors.append( "ERROR : LAST executable " + options.LASTExecutable + " does not exist ")

    if not path.exists(options.FragGeneScanExecutable):
       errors.append( "ERROR : FragGeneScan executable " + options.FragGeneScanExecutable + " does not exist ")

    if not path.exists(options.outputfolder):
       try:
          makedirs(options.outputfolder)
          makedirs(options.outputfolder + "/stage1")
          makedirs(options.outputfolder + "/stage2")
          makedirs(options.outputfolder + "/stage3")
          makedirs(options.outputfolder + "/stage4")
       except:
          errors.append( "ERROR : Output folder " + options.outputfolder + " could not be created ")

    if errors:
       for error in errors:
          print error
       return False
    return True;



def preProcess(options):
    counter = 0
    extended_file_name = options.outputfolder +"/stage1/" + options.sample_name + ".extendedFrags.fastq"
    notCombined_file_name = options.outputfolder +"/stage1/" + options.sample_name + ".notCombined.fastq"

    try:
       map_file_path = options.outputfolder +"/stage2/" + options.sample_name + ".map"
       preprocessed_file = options.outputfolder +"/stage2/" + options.sample_name + ".sofa"
       preprocessedfile = open(preprocessed_file,'w')
       mapfile = open(map_file_path, 'w')

       preProcessFile(options, extended_file_name, preprocessedfile, mapfile,  flipMe=True)

       preProcessFile(options, notCombined_file_name, preprocessedfile, mapfile)

       preprocessedfile.close()
       mapfile.close()

    except:
       print "Something went wrong"
       print traceback.print_exc(10)
       sys.exit(0)


def preProcessFile(options, filename, preprocessedfile, mapfile, flipMe=False):

    try:
       file = open(filename,'r')
    except:
       print "Error ", notCombined_file_name
       sys.exit(0) 


    regExs = [ 
               re.compile(r'@(\S+)[/]([12])'), 
               re.compile(r'@(\S+)[_]([12])'), 
               re.compile(r'@(\S+)\s([12]):[YN]:\d+:[ATCG]*')
             ]

    prev = "XXX"
    counter = 1
    seqNo  = 0
    seqName = None
    name =''
    parity=''
    for line in file:
        if counter % 4 == 1:
            success = False
            # handle regex of the 1st kind
            if success ==False:
               result = regExs[0].search(line)
               if result:
                   name = result.group(1)
                   parity = result.group(2)
                   success = True

            # handle regex of the 2nd kind
            if success==False:
               result = regExs[1].search(line)
               if result:
                   name = result.group(1)
                   parity = result.group(2)
                   success = True


            # handle regex of the 3rd kind
            if success==False:
               result = regExs[2].search(line)
               if result:
                   name = result.group(1)
                   parity = result.group(2)
                   success = True


            if name!=prev:
               seqNo +=1

            if flipMe:
               parity = '0'

            seqName = options.sample_name + "_" + str(seqNo) + "_" + parity
            prev = name
            fprintf(preprocessedfile, ">%s\n", seqName)


        if counter % 4 == 2:
            fprintf(mapfile, "%s\t%s\t%d\n", seqName, name, len(line.strip() ))
            fprintf(preprocessedfile, "%s\n", line.strip())

        counter += 1
	
    file.close()
        

def runFragGeneScan(options):
	
	
    fragGeneScan = options.FragGeneScanExecutable
    modelFile =  "illumina_10"

    inputFile =  options.outputfolder + "/stage2/" + options.sample_name + ".sofa"
    outputfile = options.outputfolder + "/stage3/" + options.sample_name 
    thread = options.threads

    arguments =  [ fragGeneScan, "-s", inputFile, "-o", outputfile, "-w", "0", "-t", modelFile, "-p", thread]
    #print ' '.join(arguments)

    result = getstatusoutput(' '.join(arguments))

    return (0, '')


def  isDuplicate(hits):
    for key in hits['1'].keys():
       if key in hits['2']:
          return True

    return False


def deDeuplicate(options, last_output_file):

    try:
        lastout = open(last_output_file, 'r')
    except:
        return;

    hits = { }

    try:
       status_file_path = options.outputfolder + PATHDELIM + "stage4" +PATHDELIM + options.sample_name + ".status"

       sofa_final_file_path = options.outputfolder + PATHDELIM + "stage4" +PATHDELIM + options.sample_name + ".SOFA.final.faa"
       no_hits_file_path = options.outputfolder + PATHDELIM + "stage4" +PATHDELIM + options.sample_name + ".nohit.faa"
       remove_dups_file_path = options.outputfolder + PATHDELIM + "stage4" +PATHDELIM + options.sample_name + ".dups.removed.faa"

       sofa_final_file = open(sofa_final_file_path, 'w')
       no_hits_file = open(no_hits_file_path, 'w')
       remove_dups_file = open(remove_dups_file_path, 'w')
      
       statusfile = open(status_file_path,'w' )
    except:
       print "Cannot open status file"
       sys.exit(0)


    sofa_final = {}
    no_hits = {}
    removed_seqs = {}


    #pattern = re.compile(r'(\S+)_([012])')
    pattern = re.compile(r'([^_\t\s]+_\d+)_([012])')
    prevContig = '---'
    for line in lastout:
         fields = [ x.strip() for x in line.split('\t') ] 
         result = pattern.search(fields[0])
         if result:
            contig = result.group(1)
            parity = result.group(2)
        
            if prevContig!=contig:
               if '1' in hits and '2' in  hits:
                   duplicate =  isDuplicate(hits)
                   fprintf(statusfile,contig +  "\t" + str(duplicate) +"\n")
                   if duplicate :
                      sofa_final[contig+"_1"] = True
                      removed_seqs[contig+"_2"] = True
                   else:
                      sofa_final[contig+"_1"] = True
                      sofa_final[contig+"_2"] = True

               elif '1' in hits and not '2' in  hits:
                   fprintf(statusfile, contig +  "\t" + 'hit1' +"\n")
                   sofa_final[contig+"_1"] = True
                   no_hits[contig+"_2"] = True
               elif not '1' in hits and  '2' in  hits:
                   fprintf(statusfile, contig +  "\t" + 'hit2' +"\n")
                   sofa_final[contig+"_2"] = True
                   no_hits[contig+"_1"] = True
               else:
                   fprintf(statusfile, contig +  "\t" + 'nohit' +"\n")
                   no_hits[contig+"_1"] = True
                   no_hits[contig+"_2"] = True

               hits = {}

            if parity =='1':
               if not '1' in hits:
                  hits['1'] = {}
               hits['1'][fields[1]] = ''

            if parity =='2':
               if not '2' in hits:
                 hits['2'] = {}
               hits['2'][fields[1]] = ''

            prevContig=contig

    statusfile.close()
    lastout.close()
    print 'Success : LAST'

 # start writing the fasta file

    namepattern = re.compile(r'([^_\t\s]+_\d+_[012])')
    fasta_file = options.outputfolder + "/stage3/" + options.sample_name + ".faa"
    fasta_reader = FastaReader(fasta_file)
    for record in fasta_reader:
        name = record.longname[1:]
        line = record.sequence

        result = namepattern.search(name)
        if result:
            extracted_name = result.group(1)
            if extracted_name in sofa_final:
              fprintf(sofa_final_file,"%s\n%s\n",record.longname, record.sequence)
            elif extracted_name in no_hits:
              fprintf(no_hits_file,"%s\n%s\n",record.longname, record.sequence)
            elif extracted_name in removed_seqs:
              fprintf(remove_dups_file,"%s\n%s\n",record.longname, record.sequence)
            else:
              fprintf(sofa_final_file,"%s\n%s\n",record.longname, record.sequence)

    sofa_final_file.close()
    no_hits_file.close()
    remove_dups_file.close()



def runDereplication(options):
    lastargs = {}
    lastargs['bitScore'] =options.bitScore
    lastargs['evalue'] =options.evalue
    lastargs['last_executable'] =options.LASTExecutable
    lastargs['refdb'] =options.refdb
    lastargs['query_file'] =  options.outputfolder + PATHDELIM + 'stage3' + PATHDELIM + options.sample_name + ".faa"
    lastargs['last_output']=  options.outputfolder + PATHDELIM + 'stage4' + PATHDELIM + options.sample_name + ".LASTout"
    _execute_LAST(lastargs )

    input = lastargs['last_output'] 
    output = lastargs['last_output'] + ".tmp"
    batch_sort(input, output, key =None, buffer_size = 5000000, tempdirs = options.tempdirs)

    deDeuplicate(options, output)

    return (0, '')

def getFragGeneScanFiles( options):

    fileMerged = options.outputfolder + PATHDELIM +  options.sample_name + ".extendedFrags.fasta"
    fileNotmerged1 = options.outputfolder + PATHDELIM +  options.sample_name + ".notCombined_1.fasta"
    fileNotmerged2 = options.outputfolder + PATHDELIM +  options.sample_name + ".notCombined_2.fasta"

    files  = [ fileMerged, fileNotmerged1, fileNotmerged2 ]
    for file in files:
       if not  path.exists(file):
          print "ERROR : Expected file " + file + " does not exist"
          return []

    return files



def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    
    
    global parser

    options, args = parser.parse_args(argv)

    if not checkInput(options):
       return 1, "ERROR"

    if '1' in options.stage: 
        print "Running : FLASH "
        result =  runFLASH(options)
        if result[0]==0:
           print "Success : FLASH "
        else:
           return (1,'ERROR : Could not run FLASH successfully!')


    if '2' in options.stage:
    	print "Running : FASTQ to FASTA "
        preProcess(options)
        if result[0]==0:
           print "Success : FASTQ to FASTA "


    if '3' in options.stage:
       print "Running : FragGeneScan+ "
       result =  runFragGeneScan(options)
       if result[0]==0:
         print "Success : FragGeneScan+ "


    if '4' in options.stage:
        print "Running : Deduplication "
        result =  runDereplication(options)
        if result[0]==0:
           print "Success : Dereplication "
        else:
          return (1,'ERROR : Could not run Dereplication successfully!')

    return (0, '')



def runFLASH(options):

    inputfile = options.inputfolder + '/' + options.sample_name + '.fastq'

    arguments =  [ options.FlashExecutable, "-O", "-I",  "-o", options.outputfolder +"/stage1/" + options.sample_name, options.inputfile ]

    #print options.FlashExecutable
    #print ' '.join(arguments)
    result = getstatusoutput(' '.join(arguments))

    return result


# this is the portion of the code that fixes the name

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes


     
def write_new_file(lines, output_file):
    
    print "Fixing file " + output_file 
    try:
       outputfile = open(output_file,'w')
       pass
    except IOError:
         print "ERROR :Cannot open output file "  + output_file
   
    for line in lines:
       fprintf(outputfile, "%s\n", line)

    outputfile.close()



import os
from tempfile import gettempdir
from itertools import islice, cycle
from collections import namedtuple
import heapq

Keyed = namedtuple("Keyed", ["key", "obj"])

def merge(key=None, *iterables):
    # based on code posted by Scott David Daniels in c.l.p.
    # http://groups.google.com/group/comp.lang.python/msg/484f01f1ea3c832d

    if key is None:
        for element in heapq.merge(*iterables):
            yield element
    else:
        keyed_iterables = [(Keyed(key(obj), obj) for obj in iterable) for iterable in iterables]
        for element in heapq.merge(*keyed_iterables):
            yield element.obj        


def batch_sort(input, output, key=None, buffer_size=32000, tempdirs=[]):

    if not tempdirs:
        tempdirs.append(gettempdir())

    chunks = []
    try:
        with open(input,'rb',64*1024) as input_file:
            input_iterator = iter(input_file)
            for tempdir in cycle(tempdirs):
                current_chunk = list(islice(input_iterator,buffer_size))
                if not current_chunk:
                    break
                current_chunk.sort(key=key)
                output_chunk = open(os.path.join(tempdir,'%06i'%len(chunks)),'w+b',64*1024)
                chunks.append(output_chunk)
                output_chunk.writelines(current_chunk)
                output_chunk.flush()
                output_chunk.seek(0)

        with open(output,'wb',64*1024) as output_file:
            output_file.writelines(merge(key, *chunks))

    finally:
        for chunk in chunks:
            try:
                chunk.close()
                os.remove(chunk.name)
            except Exception:
                pass



def MetaPathways_SOFA(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tBUILD_PGDB\n")
    createParser()
    status, error = main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (status,'')

if __name__ == '__main__':
    createParser()
    status, error = main(sys.argv[1:])

