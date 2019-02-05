#!/usr/bin/env python
"""
Make random substitutions in fasta sequences
"""

from sys    import argv,stdin,stderr,exit
from random import seed as random_seed, \
                   random as unit_random, \
                   choice as random_choice
from random_fasta_names import RandomNamer

def usage(s=None):
	message = """
usage: cat xxx.fa | random_fasta_subs <rate> [options]
  <rate>             (required) substitution rate, a probability;  this can be
                     specified as a number, percentage, or fraction
  --name=<template>  template for sequence names, as per random_fasta_names.py
  --wrap=<length>    number of characters per line
  --lower            show substitutions in lowercase
  --keepaslower      lowercase bases in input remain lowercase in output
  --seed=<string>    specify random number generator seed
  --outputboth       output the original sequences as well as the modified ones
  {<string>}         set the alphabet
                     (default is {ACGT})"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global wrapLength

	# parse the command line

	subRate      = None
	nameTemplate = None
	wrapLength   = 100
	showAsLower  = False
	keepAsLower  = False
	outputBoth   = False
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
		elif (arg == "--keepaslower"):
			keepAsLower = True
		elif (arg.startswith("--seed=")):
			random_seed(argVal)
		elif (arg == "--outputboth"):
			outputBoth = True
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (arg.startswith("{")) and (arg.endswith("}")):
			alphabet = arg[1:-1]
			if (alphabet == ""):
				usage("what am I supposed to do with an empty alphabet?")
		elif (subRate == None):
			try:
				subRate = parse_probability(arg,strict=True)
			except ValueError:
				assert (False), "\"%s\" isn't a valid probability" % arg
		else:
			usage("unrecognized option: %s" % arg)

	if (subRate == None):
		usage("you must tell me the substitution rate")

	# create the substitution table

	nucToChoices = {}
	nucToChoices[""] = [ch for ch in alphabet]

	for nuc in "ACGT":
		if (nuc in nucToChoices): continue
		nucToChoices[nuc] = [ch for ch in alphabet if (ch != nuc)]

	if (keepAsLower):
		for nuc in "ACGT".lower():
			if (nuc in nucToChoices): continue
			nucToChoices[nuc] = [ch for ch in alphabet.lower() if (ch != nuc)]

	if (showAsLower):
		for nuc in nucToChoices:
			nucToChoices[nuc] = [ch.lower() for ch in nucToChoices[nuc]]

	# process the sequences

	namer = None
	if (nameTemplate != None):
		namer = RandomNamer(nameTemplate)

	for (name,seq) in fasta_sequences(stdin):
		if (outputBoth):
			print_sequence(name,seq)

		seqLen = len(seq)
		subs = [ix for ix in xrange(seqLen) if (unit_random() < subRate)]
		if (subs != []):
			seq = [nuc for nuc in seq]
			for ix in subs:
				nuc = seq[ix]
				if (nuc in nucToChoices): choices = nucToChoices[nuc]
				else:                     choices = nucToChoices[""]
				if (choices != []): sub = random_choice(choices)
				else:               sub = random_choice(alphabet)
				seq[ix] = sub
			seq = "".join(seq)

		if (namer != None):
			name = namer.next()

		print_sequence(name,seq)


def print_sequence(name,seq):
	print ">%s" % name
	if (wrapLength == None):
		print seq
	else:
		for i in range(0,len(seq),wrapLength):
			print seq[i:i+wrapLength]


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


# parse_probability--
#	Parse a string as a probability

def parse_probability(s,strict=True):
	scale = 1.0
	if (s.endswith("%")):
		scale = 0.01
		s = s[:-1]

	try:
		p = float(s)
	except:
		try:
			(numer,denom) = s.split("/",1)
			p = float(numer)/float(denom)
		except:
			raise ValueError

	p *= scale

	if (strict) and (not 0.0 <= p <= 1.0):
		raise ValueError

	return p


if __name__ == "__main__": main()
