#!/usr/bin/env python
"""
Convert a fastq file to a file with each sequence on a single line (as four
fields plus, perhaps, a separator after the header)
"""

import sys

def main():

	# parse the command line

	separator      = "#"
	toUpper        = False
	noAt           = False
	nameParse      = None
	copyMiddleName = True

	for arg in sys.argv[1:]:
		if (arg == "--noseparator") or (arg == "--nosep"):
			separator = None
		elif (arg.startswith("--separator=")):
			separator = arg.split("=",1)[1]
		elif (arg == "--upper"):
			toUpper = True
		elif (arg == "--noat"):
			noAt = True
		elif (arg == "--nameparse=darkspace"):
			nameParse      = "darkspace"
			copyMiddleName = False
		elif (arg == "--nomiddle"):
			copyMiddleName = False
		elif (arg.startswith("--")):
			assert (False), "unknown argument: %s" % arg
		else:
			assert (False), "unknown argument: %s" % arg

	# process the sequences

	if (noAt): atSign = ""
	else:      atSign = "@"

	for (name,info) in read_fastq(sys.stdin):
		(_,seq,mid,qual) = info
		if (nameParse == "darkspace"):
			name = name.split()[0]
		if (separator != None) and (separator in name):
			name = name.replace(separator,"_")
		if (toUpper): seq = seq.upper()
		if (not copyMiddleName): mid = "+"
		if (separator == None): print "%s%s %s %s %s"    % (atSign,name,          seq,mid,qual)
		else:                   print "%s%s %s %s %s %s" % (atSign,name,separator,seq,mid,qual)


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
			   "length mismatch read vs. qualities at line %d" % lineNum
		yield (name,info)

	assert (lineNum % 4 == 0), \
		   "incomplete read at end of file"


if __name__ == "__main__": main()
