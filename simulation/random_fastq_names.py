#!/usr/bin/env python
"""
Give random names (that look like read names) to fastq sequences
"""

from sys                import argv,stdin,stderr,exit
from math               import ceil
from random             import seed as random_seed,choice
from random_fasta_names import RandomNamer


def usage(s=None):
	message = """
usage: cat xxx.fastq | random_fastq_names <name_template> [options]
  --seed=<string>            specify random number generator seed
  --head=<number>            limit the number of sequences
  --keepnames[=<separator>]  keep original names on the header line, as
                             .. <random><separator><original>
  --reportmapping            report the name mapping to stderr
  --mapping=<file>           report the name mapping to a file
  --unmapping=<file>         report the name unmapping to a file

  <template> is like this: BASE{4}_[6] where {4} is replaced by four random
  letters/numbers and [6] is replaced by the read number.

  {4}        means 4 random characters for each read.
  {4L}       means 4 random letters for each read.
  {4D}       means 4 random digits for each read.
  {{4}}      means the same 4 random characters for each read.
  {{4L}}     means the same 4 random letters for each read.
  {{4D}}     means the same 4 random digits for each read.
  [4]        meand a 4 digit counter
  [4@100001] means the counter would begin at 1001."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	template      = None
	headLimit     = None
	keepNames     = None
	reportMapping = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg == "--keepnames"):
			keepNames = " "
		elif (arg.startswith("--keepnames=")):
			assert (keepNames == None)
			keepNames = argVal
		elif (arg == "--reportmapping"):
			assert (reportMapping == None)
			reportMapping = True
		elif (arg.startswith("--unmapping=")):
			assert (reportMapping == None)
			reportMapping   = argVal
			reportAsReverse = True
		elif (arg.startswith("--mapping=")):
			assert (reportMapping == None)
			reportMapping   = argVal
			reportAsReverse = False
		elif (arg.startswith("--seed=")):
			random_seed(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (template == None):
			template = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (template == None):
		usage("you must provide a template")

	if (reportMapping == True):
		reportMappingF = stderr
	elif (reportMapping != None):
		reportMappingF = file(reportMapping,"wt")
	else:
		reportMappingF = None

	# process the sequences

	namer = RandomNamer(template)

	seqNum = 0
	for (name,info) in read_fastq(stdin):
		seqNum += 1
		if (headLimit != None) and (seqNum > headLimit):
			print >>stderr, "limit of %d sequences reached" % headLimit
			break

		name = namer.next()
		if (keepNames != None): print "@%s%s%s" % (name,keepNames,seqName)
		else:                   print "@%s"    %  name
		print "\n".join(info[1:4])

		if (reportMappingF != None):
			if (reportAsReverse):
				print >>reportMappingF, "%s %s" % (name,seqName)
			else:
				print >>reportMappingF, "%s %s" % (seqName,name)

	if (reportMappingF != None) and (reportMapping != True):
		reportMappingF.close()


# read_fastq--
#	Yield the next read from a fastq file

def read_fastq(f):
	lineNum = 0
	for line in f:
		lineNum += 1
		line = line.strip()

		if (lineNum % 4 == 1):
			assert (line.startswith("@")), \
				   "bad read name at line %d" % lineNum
			name = line[1:]
			info = [line]
			continue

		if (lineNum % 4 == 2):
			read = line
			info += [line]
			continue

		if (lineNum % 4 == 3):
			assert (line.startswith("+")), \
				   "can't understand line %d:\n%s" % (lineNum,line)
			info += [line]
			continue

		quals = line
		info += [line]
		assert (len(quals) == len(read)), \
			   "length mismatch read vs. qualities at line %d\nname: \"%s\"\nread: \"%s\"\nqual: \"%s\"\n" \
			 % (lineNum,name,read,quals)
		yield (name,info)

	if (lineNum % 4 != 0):
		message =  ["incomplete read at end of file (read %d lines)" % lineNum]
		message += ["name: \"%s\"" % name]
		if (lineNum % 4 >=2): message += ["nread: \"%s\"" % read]
		assert (False), "\n".join(message)


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()
