#!/usr/bin/env python
"""
Create a cluster job file to map reads.
"""

from sys import argv,stdin,stdout,stderr,exit


def usage(s=None):
	message = """

usage: create_script_map [options] > map.sh
  <sub>_<samp>_<type>        (required) run descriptor; for example, CS_NORM_PE
                             means subject "CS", sample "NORM", and type "PE";
                             other filenames can use "{run}" to refer to this
                             string
  --base=<path>              path prefix; other filenames can use "{base}" to
                             refer to this path
  --ref=<filename>           (required) bwa reference file
  --reads=<filename>         (required) fastq file(s);  usually this looks like
                             this:
                               {base}/reads/{run}.{mate}.fastq
  --bam=<filename>           (required) place to write bam file(s);  usually
                             this looks like this:
                               {base}/alignments/{run}
  --namesorted               create a name sorted bam file too
  --qualityfiltered          create a quality filtered bam file too
  --temp=<filename>          path to place to store temporary files
                             (default is {base}/temp/{run}.temp)
  --initialize=<text>        (cumulative) shell command to add to job beginning
                             "shebang:bash" is mapped "#!/usr/bin/env bash"
                             other commands are copied "as is" """

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global basePath,runName
	global debug

	bashShebang = "#!/usr/bin/env bash"

	# parse args

	runName            = None
	basePath           = None
	refFilename        = None
	readsFilename      = None
	bamFilename        = None
	tempFilename       = None
	doNameSorted       = False
	doQualityFiltering = False
	bashInitializers   = ["set -eu"]
	debug              = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1].strip()

		if (arg.startswith("--base=")) or (arg.startswith("--basepath=")) or (arg.startswith("--path=")):
			basePath = argVal
		elif (arg.startswith("--reference=")) or (arg.startswith("--ref=")):
			refFilename = argVal
		elif (arg.startswith("--reads=")):
			readsFilename = argVal
		elif (arg.startswith("--bam=")):
			bamFilename = argVal
		elif (arg == "--namesorted"):
			doNameSorted = True
		elif (arg == "--qualityfiltered"):
			doQualityFiltering = True
		elif (arg.startswith("--temp=")):
			tempFilename = argVal
		elif (arg.startswith("--initialize=")) or (arg.startswith("--init=")):
			if (argVal == "shebang:bash"):
				argVal = bashShebang
			if (argVal == "set -eu"):
				bashInitializers = [x for x in bashInitializers if (x != "set -eu")]
			bashInitializers += [argVal]
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (runName == None):
			fields = arg.split(":",2)
			if (len(fields) != 3):
				fields = arg.split("_")
			if (len(fields) < 3) or (fields[-1] not in ["PE","MP"]):
				usage("\"%s\" is not a valid run descriptor" % arg)
			runName = "_".join(fields)
		else:
			usage("unrecognized option: %s" % arg)

	if (runName == None):
		usage("you have to give me a run descriptor")

	if (refFilename == None):
		usage("you have to give me a bwa sequence reference")

	if (bamFilename == None):
		usage("you have to give me a bam filename")

	if (tempFilename == None):
		tempFilename = "{base}/temp/{run}.temp"

	##########
	# perform filename substitution
	##########

	if   (basePath == None):       basePath = "."
	elif (basePath.endswith("/")): basePath = basePath[:-1]

	refFilename   = do_filename_substitutition(refFilename)
	readsFilename = do_filename_substitutition(readsFilename)
	bamFilename   = do_filename_substitutition(bamFilename)
	tempFilename  = do_filename_substitutition(tempFilename)

	if (bamFilename.endswith(".bam")): bamFilename = bamFilename[:-4]

	##########
	# create the job's shell script
	##########

	# write bash intitializers

	if (bashInitializers != None):
		for (ix,bashInitializer) in enumerate(bashInitializers):
			if (bashInitializer != bashShebang): continue
			print bashInitializer
			bashInitializers[ix] = None
		for (ix,bashInitializer) in enumerate(bashInitializers):
			if (bashInitializer != "set -eu"): continue
			print bashInitializer
			bashInitializers[ix] = None
		for bashInitializer in bashInitializers:
			if (bashInitializer == None): continue
			print do_filename_substitutition(bashInitializer)
		print

	# write commands describing the files the script will create

	print "echo \"will write alignments to %s\"" % (bamFilename + ".bam")

	if (doNameSorted):
		print "echo \"will write name-sorted alignments to %s\"" % (bamFilename + ".name_sorted.bam")

	if (doQualityFiltering):
		qfBamInFilename  = bamFilename + ".bam"
		qfBamOutFilename = "%s.ql_filtered.bam" % bamFilename
		if (doNameSorted):
			qfBamInFilename  = "%s.name_sorted.bam" % bamFilename
			qfBamOutFilename = "%s.ql_filtered.name_sorted.bam" % bamFilename
		print "echo \"will write quality-filtered alignments to %s\"" % qfBamOutFilename

	# write command(s) to map the reads

	commands =  []

	command  =  ["time bwa mem"]
	command  += [refFilename]
	command  += [readsFilename.replace("{mate}","1")]
	command  += [readsFilename.replace("{mate}","2")]
	commands += [command]

	command  =  ["> " + tempFilename + ".sam"]
	commands += [command]

	print
	print "echo \"=== mapping reads ===\""
	print
	print commands_to_pipeline(commands)

	# write command(s) to position-sort the alignments and convert to bam

	commands =  []
	command  =  ["time samtools sort"]
	command  += ["-T %s" % tempFilename]
	command  += ["-o %s.bam" % bamFilename]
	command  += [tempFilename + ".sam"]
	commands += [command]

	print
	print "echo \"=== position-sorting alignments ===\""
	print
	print commands_to_pipeline(commands)

	# write command(s) to name-sort the alignments

	if (doNameSorted):
		commands =  []

		command  =  ["time samtools sort -n"]
		command  += ["-T %s" % tempFilename]
		command  += ["-o %s.name_sorted.bam" % bamFilename]
		command  += [bamFilename + ".bam"]
		commands += [command]

		print
		print "echo \"=== name-sorting alignments ===\""
		print
		print commands_to_pipeline(commands)

	# write command(s) to quality-filter the alignments

	if (doQualityFiltering):
		commands =  []

		command  =  ["time samtools view -h alignments/%s" % qfBamInFilename]
		commands += [command]

		command   = ["filtered_sam_to_intervals"]
		if (doNameSorted): command += ["--namesorted"]
		command  += ["--prohibit:\"(CIGAR == *)\""]
		command  += ["--require:\" (RNEXT == =)\""]
		command  += ["--require:\" (MAPCLIP >= RLEN*\${clipThreshold})\""]
		command  += ["--require:\" (MINCLIP <= 5)\""]
		command  += ["--require:\" (MAPQ >= \${mapQThreshold})\""]
		command  += ["--justsam"]
		command  += ["--progress=2M --progress=output:2M"]
		commands += [command]

		command  =  ["samtools view -bS -"]
		commands += [command]

		command  =  ["> alignments/%s" % qfBamOutFilename]
		commands += [command]

		print
		print "echo \"=== quality-filtering alignments ===\""
		print
		print "mapQThreshold=40"
		print "clipThreshold=0.40"
		print commands_to_pipeline(commands)

	# write command(s) to clean up

	commands =  []
	command  =  ["rm " + tempFilename + ".sam"]
	commands += [command]

	print
	print "echo \"=== cleaning up ===\""
	print
	print commands_to_pipeline(commands)


def commands_to_pipeline(commands):
	pipeline = []
	for (cmdNum,cmd) in enumerate(commands):
		if (cmdNum == 0): prefix = ""
		else:             prefix = "  | "

		if (cmd[0].startswith(">")):
			assert (cmdNum != 0)
			assert (len(cmd) == 1)
			prefix = "  "

		pipeline += [prefix + cmd[0]] 
		for line in cmd[1:]:
			pipeline += ["      " + line]

	return " \\\n".join(pipeline)


def do_filename_substitutition(s):
	if ("{base}" in s):
		assert (basePath != None)
		s = s.replace("{base}",basePath)
	if ("{run}" in s):
		assert (runName != None)
		s = s.replace("{run}",runName)
	return s


if __name__ == "__main__": main()

