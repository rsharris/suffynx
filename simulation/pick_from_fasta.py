#!/usr/bin/env python
"""
Pick a *known* subset of the sequences in a FASTA file.
"""

from sys import argv,stdin,stderr,exit

def usage(s=None):
	message = """
usage: cat fasta_file | pick_from_fasta <sequence_names_file> [options]
  <sequence_names_file>  file listing the names of sequences of interest, one-
                         per-line
  --sequence[s]=<names>  comma-separated list of sequence names
  --sort=length          output sequences in order of decreasing length
  --sort=name            output sequences in name order
  --checkduplicates      make sure that the same name is not used for different
                         sequences
  --complement           output the sequences that are NOT listed
                         (note that this is the mathematical complement of a
                         set, NOT reverse complement of a nucleotide sequence)
  --nowarn               don't warn the user when she has listed sequences that
                         are not present
  --delimiter=<char>     character that terminates the name in a fasta header
  --comma                name in a fasta header is terminated by a comma
  --trimheader           when outputing sequences, remove an excess stuff
                         from the header lines
  --notrimcomments       don't remove comments from the sequence names file;
                         they begin with a #;  note that certain sequencing
                         platforms have names with # in them
  --progress=<number>    periodically report how many lines we've read"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

    # parse the command line

	seqNamesFilename = None
	seqNamesGiven    = []
	sortBy           = None
	checkDuplicates  = False
	complement       = False
	warn             = True
	headerDelimiter  = None
	trimHeader       = False
	trimComments     = True
	reportProgress   = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--sequence=")) or (arg.startswith("--sequences=")):
			newNames = []
			for name in argVal.split(","):
				assert (name != ""), "empty name not allowed"
				if (name not in seqNamesGiven) and (name not in newNames):
					newNames += [name]
			seqNamesGiven += newNames
		elif (arg == "--sort=length"):
			sortBy = "length"
		elif (arg == "--sort=name"):
			sortBy = "name"
		elif (arg in ["--checkduplicates","--checkdups"]):
			checkDuplicates = True
		elif (arg == "--complement"):
			complement = True
		elif (arg == "--nowarn"):
			warn = False
		elif (arg.startswith("--delimiter=")):
			headerDelimiter = argVal
		elif (arg == "--comma"):
			headerDelimiter = ","
		elif (arg == "--trimheader"):
			trimHeader = True
		elif (arg == "--notrimcomments"):
			trimComments = False
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (seqNamesFilename == None):
			seqNamesFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (seqNamesFilename == None) and (seqNamesGiven == []):
		usage("you have to give me a sequence names file")

	# read the sequence names

	seqNames = {}

	for name in seqNamesGiven:
		seqNames[name] = True

	if (seqNamesFilename != None):
		lineNumber = 0
		for line in file(seqNamesFilename):
			if (trimComments): line = line.split("#")[0].strip()
			else:              line = line.strip()
			lineNumber += 1
			if (line == ""): continue
			seqNames[line] = True

	# process the fasta file, printing each desired sequence or collecting it
	# for later sort

	seqSeen     = {}
	seqToLines  = {}
	seqToHeader = {}
	currentSeq  = None

	lineNumber = seqNumber = 0
	for line in stdin:
		line = line.strip()
		lineNumber += 1

		if (line.startswith(">")):
			seqNumber += 1
			if (reportProgress != None) and (seqNumber % reportProgress == 0):
				print >>stderr, "reading sequence %s" % commatize(str(seqNumber))

			if (currentSeq != None):
				if (not checkDuplicates):
					if (sortBy == None):
						print currentHeader
						print "\n".join(currentSeq)
					else:
						seqToLines [name] = currentSeq
						seqToHeader[name] = currentHeader
				elif (name in seqToLines):
					assert ("".join(currentSeq) == "".join(seqToLines[name])), \
						   "%s is used for different sequences" % name
				else:
					seqToLines[name] = currentSeq
					if (sortBy == None):
						print currentHeader
						print "\n".join(currentSeq)
					else:
						seqToHeader[name] = currentHeader
			if (line == ">"):
				name = ""
			elif (headerDelimiter != None):
				name = line[1:].split(headerDelimiter)[0]
			else:
				name = line[1:].split()[0]
			if (name in seqNames):
				seqSeen[name] = True
			if ((name in seqNames) == complement):
				currentSeq = None
			else:
				if (trimHeader): currentHeader = ">%s" % name
				else:            currentHeader = line
				currentSeq = []
		elif (lineNumber == 1):
			assert (False), "input does not start with a fasta header line"
		elif (currentSeq != None):
			currentSeq.append(line)

	if (currentSeq != None):
		if (not checkDuplicates):
			if (sortBy == None):
				print currentHeader
				print "\n".join(currentSeq)
			else:
				seqToLines [name] = currentSeq
				seqToHeader[name] = currentHeader
		elif (name in seqToLines):
			assert ("".join(currentSeq) == "".join(seqToLines[name])), \
				   "%s is used for different sequences" % name
		else:
			seqToLines[name] = currentSeq
			if (sortBy == None):
				print currentHeader
				print "\n".join(currentSeq)
			else:
				seqToHeader[name] = currentHeader

	# if we're to sort by length, do so

	if (sortBy == "length"):
		names = []
		for name in seqToLines:
			length = sum([len(line) for line in seqToLines[name]])
			names += [(length,name)]
		names.sort()
		names.reverse()
		names = [name for (length,name) in names]

	# if we're to sort by name, do so

	elif (sortBy == "name"):
		names = [name for name in seqToLines]
		names.sort()

	# if we've sorted, output the sequences now

	if (sortBy != None):
		for name in names:
			print seqToHeader[name]
			print "\n".join(seqToLines[name])

	# warn the user about sequences that did not appear in the input

	if (warn):
		for name in seqNames:
			if (name not in seqSeen):
				print >>stderr,"warning, %s was not in input" % name


# int_with_unit--
#	Parse a strings as an integer, allowing unit suffixes

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


# commatize--
#	Convert a numeric string into one with commas.

def commatize(s):
	if (type(s) != str): s = str(s)
	(prefix,val,suffix) = ("",s,"")
	if (val.startswith("-")): (prefix,val) = ("-",val[1:])
	if ("." in val):
		(val,suffix) = val.split(".",1)
		suffix = "." + suffix

	try:    int(val)
	except: return s

	digits = len(val)
	if (digits > 3):
		leader = digits % 3
		chunks = []
		if (leader != 0):
			chunks += [val[:leader]]
		chunks += [val[ix:ix+3] for ix in xrange(leader,digits,3)]
		val = ",".join(chunks)

	return prefix + val + suffix


if __name__ == "__main__": main()
