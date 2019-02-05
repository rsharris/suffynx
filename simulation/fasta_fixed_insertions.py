#!/usr/bin/env python
"""
Make insertions of random sequences at fixed locations in fasta sequences
"""

from sys    import argv,stdin,stderr,exit
from math   import floor,ceil
from random import seed as random_seed,choice as random_choice
from random_fasta_names import RandomNamer

def usage(s=None):
	message = """
usage: cat xxx.fa | fasta_fixed_insertions <insertion_spec_file> [options]
  <insertion_spec_file>  file indicating where insertions are to be made
  --name=<template>      template for sequence names, as per random_fasta_names.py
  --wrap=<length>        number of characters per line
  --lower                show insertions in lowercase
  --seed=<string>        specify random number generator seed
  {<string>}             set the alphabet
                         (default is {ACGT})"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	specFilename = None
	nameTemplate = None
	wrapLength   = 100
	showAsLower  = False
	alphabet     = "ACGT"

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--name=")) or (arg.startswith("--names=")):
			nameTemplate = argVal
		elif (arg.startswith("--wrap=")):
			wrapLength = int(argVal)
			if (wrapLength <= 0): wrapLength = None
		elif (arg == "--lower"):
			showAsLower = True
		elif (arg.startswith("--seed=")):
			random_seed(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (arg.startswith("{")) and (arg.endswith("}")):
			alphabet = arg[1:-1]
			if (alphabet == ""):
				usage("what am I supposed to do with an empty alphabet?")
		elif (specFilename == None):
			specFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (specFilename == None):
		usage("you must give me an insertion specifications file")

	if (showAsLower):
		alphabet = alphabet.lower()

	# read the insertion specifications

	nameToInsertions = {}

	f = file(specFilename,"rt")

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) >= 3), \
		       "not enough fields at line %d (%d, expected at least 3)" \
		     % (lineNumber,len(fields))

		name   = fields[0]
		start  = parse_insertion_value(fields[1])
		length = parse_insertion_value(fields[2])

		if (name not in nameToInsertions): nameToInsertions[name] = []
		nameToInsertions[name] += [(start,length)]

	for name in nameToInsertions:
		nameToInsertions[name].sort()
		nameToInsertions[name].reverse()

	f.close()

	# process the sequences

	nameSeen = {}

	namer = None
	if (nameTemplate != None):
		namer = RandomNamer(nameTemplate)

	for (name,seq) in fasta_sequences(stdin):
		if (name in nameToInsertions):
			nameSeen[name] = True
			for (start,length) in nameToInsertions[name]:
				seqLen = len(seq)
				if (type(start)  == float): start  = int(floor(start*seqLen))
				if (type(length) == float): length = int(floor(length*seqLen))
				seq = seq[:start] \
				    + random_sequence(length,alphabet=alphabet) \
				    + seq[start:]

		if (namer != None):
			name = namer.next()

		print ">%s" % name
		if (wrapLength == None):
			print seq
		else:
			seqLen = len(seq)
			for i in range(0,seqLen,wrapLength):
				print seq[i:i+wrapLength]

	for name in nameToInsertions:
		if (name in nameSeen): continue
		print >>stderr, "WARNING: sequence \"%s\" was not in the input" % name


# random_sequence--

def random_sequence(seqLen,alphabet=None):
	if (alphabet == None): alphabet = ["A","C","G","T"]
	elif (type(alphabet) == str): alphabet = [ch for ch in alphabet]

	return "".join([random_choice(alphabet) for ix in xrange(seqLen)])


# fasta_sequences--
#	Read the fasta sequences from a file

def fasta_sequences(f):
	seqName = None
	seqNucs = None

	for line in f:
		line = line.strip()
		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seqNucs))
			seqName = line[1:].strip()
			seqNucs = []
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += [line]

	if (seqName != None):
		yield (seqName,"".join(seqNucs))


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


# parse_insertion_value--
#	These return an int (indicating an exact position) or a float (indicating
#	a relative position

def parse_insertion_value(s):
	# try parsing it as an int first

	try: return int_with_unit(s)
	except ValueError: pass

	# since it isn't an int, try parsing it as a float

	scale = 1.0
	if (s.endswith("%")):
		scale = 0.01
		s = s[:-1]

	try: return scale * float(s)
	except: pass

	try:
		(numer,denom) = s.split("/",1)
		numer = float(numer)
		denom = float(denom)
		return (scale * numer)/denom
	except: pass

	raise ValueError


if __name__ == "__main__": main()
