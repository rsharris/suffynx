#!/usr/bin/env python
"""
Randomly rearrange the lines of a file
"""

from sys    import argv,stdin,stdout,stderr,exit
from math   import ceil
from random import seed as random_seed,randint


def usage(s=None):
	message = """

usage: cat file | shuffle [seed] [options]
  --seed=<string>            set random seed
  #<string>                  set random seed
  <seed>                     set random seed, interpreted as an integer
  --block=<number>           shuffle each continguious block of lines
                             independently (less memory but not as random)
  --limit=<number>           limit the number of shiffled lines output
  --passcomments[=<prefix>]  copy comment lines directly to output; the default
                             <prefix> is a #"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global passComments

	# parse the command line

	seed         = None
	blockSize    = None
	lineLimit    = None
	passComments = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--seed=")):
			seed = argVal
		elif (arg.startswith("#")):
			seed = arg[1:]
		elif (arg.startswith("--block=")):
			blockSize = int_with_unit(argVal)
		elif (arg.startswith("--limit=")):
			lineLimit = int_with_unit(argVal)
		elif (arg == "--passcomments"):
			passComments = ["#"]
		elif (arg.startswith("--passcomments=")):
			passComments = argVal
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			try:
				seed = int(arg)
			except:
				# this is to prevent the common usage mistake of putting a filename
				# in arg 1;  the user who names her file as an integer unfortunately
				# loses this protection
				usage ("Hey! Seed must be an integer")

	if (seed != None): random_seed(seed)

	# collect lines

	if (blockSize == None):
		data = read_block()
		write_shuffled_block(data,lineLimit=lineLimit)
	else:
		while (True):
			data = read_block(numLines=blockSize)
			numLines = len(data)
			if (numLines == 0): break

			if   (lineLimit == None):    blockLimit = None
			elif (lineLimit > numLines): blockLimit = numLines
			else:                        blockLimit = lineLimit
			write_shuffled_block(data,lineLimit=blockLimit)

			if (lineLimit != None):
				lineLimit -= blockLimit
				if (lineLimit == 0): break


def read_block(numLines=None):
	data = []

	lineNum = 0
	for line in stdin:
		lineNum += 1
		line = line.rstrip()
		if (passComments != None) and (line[0] in passComments):
			print line
			continue
		data += [line]

		if (numLines != None) and (lineNum == numLines): break

	return data


def write_shuffled_block(data,lineLimit=None):
	if (lineLimit == None): lineLimit = len(data)

	for ix in range(0,lineLimit):
		# choose a random line from those remaining
		randIx = randint(ix,len(data)-1)
		line   = data[randIx]
		# 'swap' it with this one;  for completeness we'd need data[ix] = line
		data[randIx] = data[ix]
		# print it
		print line


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
