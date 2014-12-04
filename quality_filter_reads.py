#!/usr/local/env python

##print __name__

import os
import optparse

usage_line = """
quality_filter_reads.py

Version 1.0 (3 December, 2014)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

Script that quality filters raw fastq read files using the Trimmomatic program. The script allows \
for the full functionality of Trimmomatic v. 0.32, but has appropriate defaults in place for quick \
and easy running. User must specify whether the reads are single-end ("-s") or paired-end ("-p"), \
the directory (without beginning or ending /) containing the read files (with no other file types), \
and the file extensions for the read files (usually "fastq"). Naming conventions for read files are \
as follows:
1. The file extension (normally fastq) must match that passed to the command with the '-e' flag.
2. Prior to the extension, there must be an indication of the read type (single or paired) as follows:
        a. P1 and P2 for paired reads, with P1 designated for the single-end reads and P2 designated for the paired-end reads
        b. S1 for non-paired, single-end reads
3. A file root with the name of the sample
Example: ID1234_Loc1_ACTTAG-GTACAG.P1.fastq = single-end reads of paired reads of sample ID1234_Loc1_ACTTAG-GTACAG
Note: Reads that do not have this file name formatting will probably not be run correctly.

The user can also specify the number of threads to use and various parameters and settings to be \
used by Trimmomatic for read quality filtering. User must select whether quality trimming will be done \
using the default sliding window approach ("--sliding_window") or using a maxinfo algorithm ("--maxinfo"). \ 
Default values for all parameters should be sufficient for most applications. 

Trimming steps take place in the following order (assuming all options are used):
1. Bases are removed from the start of the reads (--head_crop)
2. Bases are removed from the end of the reads (--tail_crop)
3. Adapters are trimmed based on input fasta and parameters (--adapters)
4. Bases below the designated quality score are removed from the start of the reads (--leading)
5. Bases below the designated quality score are removed from the end of the reads (--trailing)
6. Reads are quality trimmed using either a sliding window (--sliding_window) or maxinfo (--maxinfo) approach
7. Reads below a designed length are removed (--drop)

See the user manual at http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf \
for more information and guidance.

The only dependency for this script is Trimmomatic v. 0.32, which must be installed in the users path.

python quality_filter_reads.py -t <#threads> -p/-s -d <directory> -e <read_file_extension> \
[--tail_crop <length> --head_crop <length> --adapters <adapters.fasta> --adapters_mismatch <#> \
--adapters_palindrome <#> --adapters_simple <#> --adapters_length <#> --adapters_keep <true/false> \
--leading <#> --trailing <#> --sliding_window --sliding_width <#> --sliding_score <#> --maxinfo \
--maxinfo_length <#> maxinfo_strictness <0.0-1.0> --drop <length> --phred <33/64>
"""

#################################################
###           Parse command options           ###
#################################################

usage = usage_line

parser = optparse.OptionParser(usage=usage)
parser.add_option("-t", action="store", type = "string", dest = "threads", help = "threads to use [2]", default = "2")
parser.add_option("-p", action="store_true", dest = "paired", help = "paired reads flag")
parser.add_option("-s", action="store_true", dest = "single", help = "single-end reads flag")
parser.add_option("-d", action="store", type = "string", dest = "directory", help = "directory containing read files")
parser.add_option("-e", action="store", type = "string", dest = "ext", help = "the file extension of the read files")
parser.add_option("--tail_crop", action = "store", type = "string", dest = "tailcrop", help = "crop read to designated length by trimming from end of reads [off]")
parser.add_option("--head_crop", action = "store", type = "string", dest = "headcrop", help = "crop designated number of bases from beginning of reads [off]")
parser.add_option("--adapters", action = "store", type = "string", dest = "adapters", help = "fasta file with adapters to be trimmed [none]")
parser.add_option("--adapters_mismatch", action = "store", type = "string", dest = "adapt_mis", help = "mismatches allowed in adapter searching seed [2]", default = "2")
parser.add_option("--adapters_palindrome", action = "store", type = "string", dest = "adapt_palin", help = "threshold for match between two adapter ligated reads for PE palinrome read alignment [30]", default = "30")
parser.add_option("--adapters_simple", action = "store", type = "string", dest = "adapt_simple", help = "accuracy of match between adapter and read sequences [10]", default = "10")
parser.add_option("--adapters_length", action = "store", type = "string", dest = "adapt_length", help = "minimum adapter length detection threshold [8]", default = "8")
parser.add_option("--adapters_keep", action = "store", type = "string", dest = "adapt_keep", help = "keep both reads during palindrome mode detection [false]", default = "false")
parser.add_option("--leading", action = "store", type = "string", dest = "leading", help = "minimum quality score to keep leading bases [10]", default = "10")
parser.add_option("--trailing", action = "store", type = "string", dest = "trailing", help = "minimum quality score to keep trailing bases [10]", default = "10")
parser.add_option("--sliding_window", action = "store_true", dest = "slid_wind", help = "use a sliding window appraoch for quality trimming instead of the maxinfo algorithm [true]", default = True)
parser.add_option("--sliding_width", action = "store", type = "string", dest = "win_width", help = "sliding window width (bp) [4]", default = "4")
parser.add_option("--sliding_score", action = "store", type = "string", dest = "win_score", help = "average quality per base score threshold for sliding window [15]", default = "15")
parser.add_option("--maxinfo", action = "store_true", dest = "maxinfo", help = "use maxinfo algorithm for quality trimming instead of a sliding window approach [false]", default = False)
parser.add_option("--maxinfo_length", action = "store", type = "string", dest = "maxinfo_length", help = "maxinfo target length (bp) parameter for maxinfo algorithm [75]", default = "75")
parser.add_option("--maxinfo_strictness", action = "store", type = "string", dest = "maxinfo_strict", help = "maxinfo strictness parameter for maxinfo algorithm [0.8]", default = "0.8")
parser.add_option("--drop", action = "store", type = "string", dest = "drop", help = "drop reads below a certain length [36]", default = "36")
parser.add_option("--phred", action = "store", type = "string", dest = "phred", help = "PHRED encoding for output read quality scores [33]", default = "33")

options, args = parser.parse_args()


#################################################
###         Setup mapping enviornment         ###
#################################################

def setup():
        os.system("mkdir filtered")                                                                                      # make 'filtered' directory (may error if already present)


#################################################
###         Filter single-end reads           ###
#################################################

def SE_filter(name, adapters, tail, head, algorithm):
	adapters = adapters
	tail = tail
	head = head
	algorithm = algorithm
	foo = name.split(".")                                                                                    # split by '.'
	print "\n***Quality filtering single-end reads from "+foo[0]+"***\n"
	input = foo[0]+".SE."+options.ext                                                               # input SE file for filtering
	SEclean = "trimmomatic-0.32.jar SE -threads "+options.threads+" -trimlog ./filtered/"+foo[0]+".qtrim.log ./raw/"+input+" ./filtered/"+foo[0]+".S1.qtrim "+head+tail+adapters+"LEADING:"+options.leading+" TRAILING:"+options.trailing+" "+algorithm+"MINLEN:"+options.drop+" TOPHRED"+options.phred	
	print SEclean
	os.system(str(SEclean))


#################################################
###         Filter paired-end reads           ###
#################################################

def PE_filter(name, adapters, tail, head, algorithm):
        adapters = adapters
        tail = tail
        head = head
        algorithm = algorithm
        foo = name.split(".")                                                                                    # split by '.'
        print "\n***Quality filtering single-end reads from "+foo[0]+"***\n"
        forward = foo[0]+".P1."+options.ext                                                               # input forward PE file for filtering
	reverse = foo[0]+".P2."+options.ext								  # input reverse PE file for filtering
        PEclean = "trimmomatic-0.32.jar PE -threads "+options.threads+" -trimlog ./filtered/"+foo[0]+".qtrim.log ./raw/"+forward+" ./raw/"+reverse+" ./filtered/"+foo[0]+".P1.qtrim ./filtered/"+foo[0]+".S1.qtrim ./filtered/"+foo[0]+".P2.qtrim ./filtered/"+foo[0]+".S2.qtrim "+head+tail+adapters+"LEADING:"+options.leading+" TRAILING:"+options.trailing+" "+algorithm+"MINLEN:"+options.drop+" TOPHRED"+options.phred
        print PEclean
        os.system(str(PEclean))


#################################################
###                Full Program               ###
#################################################

def main():
        if options.directory is None:                                                                   # User didn't input directory containing read files
                print "\n***Error: specify directory containing read files!***\n"
        if options.ext is None:                                                                                 # User didn't specify read file extension (e.g., fastq)
                print "\n***Error: specify the file extension for the read files!***\n"
        else:
                setup()                                                                                                         # Setup the pipeline environment
                print "\n***Running quality trimming pipeline***\n"
                for root,dirs,files in os.walk(options.directory):                      # Gather all objects in specified directory
                        print "\n***Gathing read files from specified directory***\n"
                names = {}
                for file in files:
                        if file.endswith("."+options.ext):                                              # If file ends with specified extension
                                foo = file.split(os.extsep)
                                name = foo[0]                                                                           # Take root of file name (everything up to 1st period)
                                if name not in names.keys():
                                        names[name] = 1                                                                 # Store each unique file root in dictionary
#               print names
		if options.adapters is None:
			print "\n***No adapter sequences specified for adapter trimming!***\n"
			adapters = ""
		else:
			adapters = "ILLUMINACLIP:"+options.adapters+":"+options.adapt_mis+":"+options.adapt_palin+":"+options.adapt_simple+":"+options.adapt_length+" "
		if options.tailcrop is None:
			print "\n***The end of the read will not be trimmed!***\n"
			tail = ""
		else:
			tail = "CROP:"+options.tailcrop+" "
		if options.headcrop is None:
			print "\n***The beginning of the read will not be trimmed!***\n"
			head = ""
		else:
			head = "HEADCROP:"+options.headcrop+" "
		if options.maxinfo == True:
			options.slid_wind == False
			algorithm = "MAXINFO:"+options.maxinfo_length+":"+options.maxinfo_strict+" "
		elif options.slid_wind == True:
			algorithm = "SLIDINGWINDOW:"+options.win_width+":"+options.win_score+" "
		else:
			print "\n***Error: specify either sliding window or maxinfo trimming with parameters!***\n"
                for name in names.keys():                                                                       # For each unique file name in dictionary
                        if options.paired == True:
				PE_filter(name, adapters, tail, head, algorithm)
				print "\n***Quality filtering complete! See 'filtered' directory for results!***\n"
			elif options.single == True:
				SE_filter(name)
				print "\n***Quality filtering complete! See 'filtered' directory for results!***\n"
			else:
				print "\n***Error: specify whether reads are single-end only ('-s') or paired end ('-p')!***\n"


#################################################
###              Run Full Program             ###
#################################################

main()
