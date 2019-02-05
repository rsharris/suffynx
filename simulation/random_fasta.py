#!/usr/bin/env python
"""
Create random fasta sequences
"""

from sys    import argv,stdin,stderr,exit
from math   import ceil
from random import seed as random_seed,choice as random_choice
from random_fasta_names import RandomNamer

def usage(s=None):
	message = """
usage: [cat lengths] | random_fasta [<number>x]<length> [options]
  <number>           number of sequences to generate
  <length>           (required) length of each sequence;  if this is "stdin",
                     lengths are read from input
  --name=<template>  template for sequence names, as per random_fasta_names.py
  --noname           do not give sequences names
  --wrap=<length>    number of characters per line
  --seed=<string>    specify random number generator seed
  {<string>}         set the alphabet
                     (default is {ACGT})"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	numSequences   = None
	sequenceLength = None
	nameTemplate   = None
	wrapLength     = 100
	alphabet       = "ACGT"

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--name=")) or (arg.startswith("--names=")):
			nameTemplate = argVal
		elif (arg == "--noname") or (arg == "--nonames"):
			nameTemplate = ""
		elif (arg.startswith("--wrap=")):
			wrapLength = int(argVal)
			if (wrapLength <= 0): wrapLength = None
		elif (arg.startswith("--seed=")):
			random_seed(argVal)
		elif (arg.startswith("--n=")) or (arg.startswith("N=")):
			numSequences   = int_with_unit(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (arg.startswith("{")) and (arg.endswith("}")):
			alphabet = arg[1:-1]
			if (alphabet == ""):
				usage("what am I supposed to do with an empty alphabet?")
		elif (sequenceLength == None):
			if ("x" in arg):
				(n,l) = arg.split("x",1)
				numSequences   = int_with_unit(n)
				sequenceLength = l
			else:
				sequenceLength = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (sequenceLength == None):
		usage("you must tell me the sequence length")
	elif (sequenceLength == "stdin"):
		pass
	else:
		sequenceLength = int_with_unit(sequenceLength)
		if (numSequences == None): numSequences = 1

	if (nameTemplate == None):
		if (numSequences == None):
			numDigits = 9
		else:
			numDigits = 1
			biggestNumber = 9
			while (biggestNumber < numSequences):
				numDigits += 1
				biggestNumber = 10*biggestNumber + 9
		nameTemplate = "RAND_[%d]" % numDigits

	# generate the sequences

	if (nameTemplate == ""): namer = None
	else:                    namer = RandomNamer(nameTemplate)

	lineNum = seqNum = 0
	while (True):
		lineNum += 1
		if (numSequences != None) and (seqNum >= numSequences): break

		seqLen = sequenceLength
		if (sequenceLength == "stdin"):
			line = stdin.readline()
			if (line == ""): break
			line = line.strip()
			if (line == ""): continue
			try:
				seqLen = int_with_unit(line)
				if (seqLen < 0): raise ValueError
			except ValueError:
				assert (False), "bad length (line %d): %s" % (lineNum,line)
			if (seqLen == 0): continue

		seqNum += 1
		seq = random_sequence(seqLen,alphabet=alphabet)

		if (namer != None):
			print ">%s" % namer.next()

		if (wrapLength == None):
			print seq
		else:
			for i in range(0,seqLen,wrapLength):
				print seq[i:i+wrapLength]


# random_sequence--

def random_sequence(seqLen,alphabet=None):
	if (alphabet == None): alphabet = ["A","C","G","T"]
	elif (type(alphabet) == str): alphabet = [ch for ch in alphabet]

	return "".join([random_choice(alphabet) for ix in xrange(seqLen)])


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
