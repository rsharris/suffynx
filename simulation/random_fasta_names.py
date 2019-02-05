#!/usr/bin/env python
"""
Give random names (that look like read names) to fasta sequences
"""

from sys    import argv,stdin,stderr,exit
from math   import ceil
from random import seed as random_seed,choice
from re     import compile


def usage(s=None):
	message = """
usage: cat xxx.fa | random_fasta_names <name_template> [options]
  --seed=<string>            specify random number generator seed
  --head=<number>            limit the number of sequences
  --keepnames[=<separator>]  keep original names on the header line, as
                             .. <random><separator><original>
  --reportmapping            report the name mapping to stderr
  --mapping=<file>           report the name mapping to a file
  --unmapping=<file>         report the name unmapping to a file

  <template> is like this: BASE{4}_[6] where {4} is replaced by four random
  letters/numbers and [6] is replaced by the read number.

  {4}      means 4 random characters for each read.
  {4L}     means 4 random letters for each read.
  {4D}     means 4 random digits for each read.
  {{4}}    means the same 4 random characters for each read.
  {{4L}}   means the same 4 random letters for each read.
  {{4D}}   means the same 4 random digits for each read.
  [4]      means a 4 digit counter
  [4@1001] means the counter would begin at 1001."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	template        = None
	headLimit       = None
	keepNames       = None
	reportMapping   = None
	reportAsReverse = False

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
	for (seqName,seq) in read_fasta_sequences(stdin):
		seqNum += 1
		if (headLimit != None) and (seqNum > headLimit):
			print >>stderr, "limit of %d sequences reached" % headLimit
			break

		name = namer.next()
		if (keepNames != None): print ">%s%s%s" % (name,keepNames,seqName)
		else:                   print ">%s"    %  name
		print "\n".join(seq)

		if (reportMappingF != None):
			if (reportAsReverse):
				print >>reportMappingF, "%s %s" % (name,seqName)
			else:
				print >>reportMappingF, "%s %s" % (seqName,name)

	if (reportMappingF != None) and (reportMapping != True):
		reportMappingF.close()


# RandomNamer--
#	Class to generate 'random' names according to a template.

class RandomNamer(object):

	def __init__(self,template,duplicatesOK=False,attemptsLimit=10000):
		self.template      = template
		self.attemptsLimit = attemptsLimit

		# determine which tokens are present in the template

		self.randomized = []
		for w in xrange(20,0,-1):
			brack = "{{%d}}" % w			# e.g. "{{5}}" means the same 5
			fields = template.split(brack)	# .. random characters for all reads
			if (len(fields) > 1):
				self.randomized += [(w,brack,"",self.rand_chars(w))]

			brack = "{{%dL}}" % w			# e.g. "{{5L}}" means the same 5
			fields = template.split(brack)	# .. random letters for all reads
			if (len(fields) > 1):
				self.randomized += [(w,brack,"",self.rand_letters(w))]

			brack = "{{%dD}}" % w			# e.g. "{{5D}}" means the same 5
			fields = template.split(brack)	# .. random digits for all reads
			if (len(fields) > 1):
				self.randomized += [(w,brack,"",self.rand_digits(w))]

			brack = "{%d}" % w				# e.g. "{5}" means 5 random
			fields = template.split(brack)	# .. characters, "different" for
			if (len(fields) > 1):			# .. each read
				self.randomized += [(w,brack,"R",None)]

			brack = "{%dL}" % w				# e.g. "{5L}" means 5 random
			fields = template.split(brack)	# .. letters, "different" for each
			if (len(fields) > 1):			# .. read
				self.randomized += [(w,brack,"L",None)]

			brack = "{%dD}" % w				# e.g. "{5D}" means 5 random
			fields = template.split(brack)	# .. digits, "different" for each
			if (len(fields) > 1):			# .. read
				self.randomized += [(w,brack,"D",None)]

		self.counted     = None				# e.g. "[5]" means 5 digit count
		self.countOffset = 0				# or "[5@10001] means start at 10001
		for w in xrange(13,-1,-1):			# ... but only longest is accepted
			if (w > 0):
				brack   = "[%d]" % w
				brackRe = compile("\[%d@(?P<start>(0|[1-9][0-9]*))\]" % w)
				fmt     = "%%0%dd" % w
			else:
				brack   = "[]"
				brackRe = compile("\[@(?P<start>(0|[1-9][0-9]*))\]")
				fmt     = "%d"
			fields = template.split(brack,1)
			if (len(fields) > 1):
				self.counted = (w,brack,fmt)
				break
			m = brackRe.search(template)
			if (m != None):
				start = m.group("start")
				self.countOffset = int(start)-1
				match = m.group(0)
				fields = template.split(match,1)
				self.template = fields[0] + brack + fields[1]
				self.counted = (w,brack,fmt)
				break

		# if necessary, create a hash to prevent duplicate names

		if (duplicatesOK) or (self.randomized == []):
			self.nameUsed = None
		else:
			self.nameUsed = {}

		self.count = 0


	def next(self):
		self.count += 1
		name = None
		for attempt in xrange(self.attemptsLimit):
			cand = self.generate_name()
			if (self.nameUsed == None) or (cand not in self.nameUsed):
				name = cand
				break
		assert (name != None), "failed to generate name %d for %s" \
		                     % (self.count,self.template)
		if (self.nameUsed != None):
			self.nameUsed[name] = True

		return name


	def generate_name(self):
		name = self.template

		# substitute, e.g. "{5}" with 5 random characters

		for (w,brack,kind,val) in self.randomized:
			fields = name.split(brack)
			name = fields[0]
			if (kind == ""):
				for f in fields[1:]:
					name += val + f
			elif (kind == "L"):
				for f in fields[1:]:
					name += self.rand_letters(w) + f
			elif (kind == "D"):
				for f in fields[1:]:
					name += self.rand_digits(w) + f
			else: # if (kind == "R"):
				for f in fields[1:]:
					name += self.rand_chars(w) + f

		# substitute, e.g. "[5]" with "%05d" % count

		if (self.counted != None):
			(w,brack,fmt) = self.counted
			fields = name.split(brack,1)
			name = (fmt % (self.count+self.countOffset)).join(fields)

		return name


	def rand_chars(self,n):
		return "".join([choice("ACDEFGHJKLMNPQRTUVWXYZ0123456789") 
		                                                  for i in range(n)])

	def rand_letters(self,n):
		return "".join([choice("ACDEFGHJKLMNPQRTUVWXYZ") for i in range(n)])

	def rand_digits(self,n):
		return "".join([choice("0123456789") for i in range(n)])


# read_fasta_sequences--

def read_fasta_sequences(f):
	name  = "(no_name)"
	lines = []

	lineNum = 0
	for line in f:
		lineNum += 1
		line = line.rstrip()

		if (line.startswith(">")):
			if (lines != []): yield (name,lines)
			name  = line[1:].strip()
			lines = []
		else:
			lines += [line]

	if (lines != []): yield (name,lines)


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
