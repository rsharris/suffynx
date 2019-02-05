#!/usr/bin/env python
"""
Merge several fasta sequences into one, with runs of Xs between them.
"""

from sys import argv,stdin,stderr


def usage(s=None):
	message = """
usage: cat xxx.fa | fasta_X_merge [options]
  --fasta=<name>          name for the resulting sequence's header
                          (by default we write no fasta header)
  --subsample=<k>/<n>     keep only the <k>th sequence of every group of <n>
                          sequences;  <k> ranges from 1 to <n>
  --limit=<number>        limit the number of input sequences read
  --separator=<length>    the number of Xs that will separate the input
                          sequences in the combined output sequence    
  --separator=none        don't separate the sequences
  --separator=N           use N as a separator instead of X
  --wrap=<length>         number of characters per line
  --progress=<count>      report status every <count> sequences
  --makeindex=<filename>  create an index file (see note below)
  --origin=one            index intervals are origin-one, closed
  --origin=zero           index intervals are origin-zero, half-open

Note:  When --makeindex is used, the resulting file maps sequence name to the
corresponding interval in the combined sequence."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	fastaName       = None
	subsampleK      = None
	subsampleN      = None
	sequenceLimit   = None
	separatorLength = 50
	separator       = "X"
	wrapLength      = 100
	reportProgress  = None
	indexFilename   = None
	origin          = "zero"
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--fasta=")) or (arg.startswith("--name=")):
			fastaName = argVal
		elif (arg.startswith("--subsample=")):
			val = argVal
			(k,n) = val.split("/",2)
			subsampleK = int(k)
			subsampleN = int(n)
			assert (0 < subsampleK <= subsampleN)
		elif (arg.startswith("--limit=")):
			sequenceLimit = int(argVal)
			assert (sequenceLimit >= 0)
			if (sequenceLimit == 0): sequenceLimit = None
		elif (arg == "--separator=none") or (arg == "--sep=none"):
			separatorLength = None
		elif (arg == "--separator=N") or (arg == "--sep=N") \
		  or (arg == "--separator=n") or (arg == "--sep=n"):
			separator = "N"
		elif (arg.startswith("--separator=")) or (arg.startswith("--sep=")):
			separatorLength = int(argVal)
			assert (separatorLength >= 0)
			if (separatorLength == 0): separatorLength = None
		elif (arg.startswith("--wrap=")):
			wrapLength = int(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int(argVal)
			assert (reportProgress >= 0)
			if (reportProgress == 0): reportProgress = None
		elif (arg.startswith("--makeindex=")) or (arg.startswith("--index=")):
			indexFilename = argVal
		elif (arg.startswith("--origin=")):
			origin = argVal
			if (origin == "0"): origin = "zero"
			if (origin == "1"): origin = "one"
			assert (origin in ["zero","one"]), "can't understand %s" % arg
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the sequences

	indexF = None
	if (indexFilename != None):
		indexF = file(indexFilename,"wt")

	offset = 0
	numSequences = 0
	outSequences = 0
	combinedSeq = ""
	for (name,seq) in fasta_sequences(stdin):
		numSequences += 1
		if (reportProgress != None) and (numSequences % reportProgress == 0):
			print >>stderr, "%s %d" % (name,numSequences)
		if (subsampleN != None):
			if ((numSequences-1) % subsampleN != (subsampleK-1)):
				continue
		if (len(seq) == 0): continue

		outSequences += 1
		if (sequenceLimit != None) and (outSequences > sequenceLimit):
			print >>stderr, "output sequence limit of %d reached" \
			              % sequenceLimit
			break

		if (offset == 0):
			if (fastaName != None): print ">%s" % fastaName
		elif (separatorLength == None):
			pass
		else:
			combinedSeq += separatorLength * separator
			offset += separatorLength

		if (indexF != None):
			(start,end) = (offset,offset+len(seq))
			if (origin == "one"): start += 1
			print >>indexF, "%s\t%d\t%d" % (name,start,end)

		combinedSeq += seq
		numLines    = len(combinedSeq) / wrapLength
		writeLen    = numLines*wrapLength
		writeSeq    = combinedSeq[:writeLen]
		combinedSeq = combinedSeq[writeLen:]
		for i in range(0,writeLen,wrapLength):
			print "".join(writeSeq[i:i+wrapLength])

		offset += len(seq)

	if (combinedSeq != ""):
		print combinedSeq

	if (indexF != None):
		indexF.close()


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


if __name__ == "__main__": main()
