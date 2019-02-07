#!/usr/bin/env python
"""
Create a cluster job file to create a clipped-breakpoints track.
"""

from sys import argv,stdin,stdout,stderr,exit


def usage(s=None):
	message = """

usage: create_script_clipped_breakpoints [options] > clipped_breakpoints.sh
  <sub>_<samp>_<type>        (required) run descriptor; for example, CS_NORM_PE
                             means subject "CS", sample "NORM", and type "PE";
                             other filenames can use "{run}" to refer to this
                             string
  --control=<filename>       read control values from a file (see list below)
  --base=<path>              path prefix; other filenames can use "{base}" to
                             refer to this path
  --bam=<filename>           (required) bam file to process
  --namesorted               the bam file has been sorted by read names
  --chromosomes=<filename>   read chromosome names and lengths from a file
                             (default is {base}/data/hg19.chrom_lengths)
  --track=<filename>         (required) track file to create
                             (default is {base}/tracks/{run}.clipped_breakpoints)
  --undated                  don't include today's date in the track name
  --bigwig[=<filename>]      create bigwig file in addition to track file
  --bigwigchroms=<filename>  chromosomes file for bedGraphToBigWig
                             (default is {base}/temp/ucsc.hg19.chrom_lengths)
  --bigwigurl=<url>          url for the bigwig file; this can use {bigwig}
                             for the bigwig filename
  --bigwiglink=<filename>    path at which to create a symbolic link to the
                             bigwig and info files; this can use {bigwig}
                             for the bigwig filename
  --bigwigposition=<interval> intitial UCSC browser interval for bigwig track
  --initialize=<text>        (cumulative) shell command to add to job beginning
                             "shebang:bash" is mapped "#!/usr/bin/env bash"
                             other commands are copied "as is"
  --head=<number>            limit the number of bam records

values read from control file:
  clipThreshold
  clipped_breakpoints.windowSize"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global basePath,runName
	global debug

	bashShebang = "#!/usr/bin/env bash"

	# parse args

	runName              = None
	controlFilename      = None
	basePath             = None
	bamFilename          = None
	isNameSorted         = False
	chromsFilename       = None
	trackName            = None
	dateInTrackname      = True
	bigWigFilename       = None
	bigWigChromsFilename = None
	bigWigUrl            = None
	bigWigLink           = None
	bigWigPosition       = None
	bashInitializers     = ["set -eu"]
	headLimit            = None
	debug                = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1].strip()

		if (arg.startswith("--control=")):
			controlFilename = argVal
		elif (arg.startswith("--base=")) or (arg.startswith("--basepath=")) or (arg.startswith("--path=")):
			basePath = argVal
		elif (arg.startswith("--bam=")):
			bamFilename = argVal
		elif (arg == "--namesorted"):
			isNameSorted = True
		elif (arg.startswith("--chromosomes=")) or (arg.startswith("--chroms=")):
			chromsFilename = argVal
		elif (arg.startswith("--track=")):
			trackName = argVal
		elif (arg == "--undated"):
			dateInTrackname = False
		elif (arg == "--bigwig"):
			bigWigFilename = "{track}.bw"
		elif (arg.startswith("--bigwig=")):
			bigWigFilename = argVal
		elif (arg.startswith("--bigwigchromosomes=")) or (arg.startswith("--bigwigchroms=")):
			bigWigChromsFilename = argVal
		elif (arg.startswith("--bigwigurl=")) or (arg.startswith("--url=")):
			bigWigUrl = argVal
		elif (arg.startswith("--bigwiglink=")) or (arg.startswith("--link=")):
			bigWigLink = argVal
		elif (arg.startswith("--bigwigposition=")) or (arg.startswith("--bigwigpos=")):
			bigWigPosition = argVal
		elif (arg.startswith("--initialize=")) or (arg.startswith("--init=")):
			if (argVal == "shebang:bash"):
				argVal = bashShebang
			if (argVal == "set -eu"):
				bashInitializers = [x for x in bashInitializers if (x != "set -eu")]
			bashInitializers += [argVal]
		elif (arg.startswith("--head=")):
			headLimit = argVal
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

	if (controlFilename == None):
		usage("you have to give me a control filename")

	if (bamFilename == None):
		usage("you have to give me a bam filename")

	if (chromsFilename == None):
		chromsFilename = "{base}/data/hg19.chrom_lengths"

	if (trackName == None):
		trackName = "{base}/tracks/{run}.clipped_breakpoints"

	if (bigWigFilename != None):
		if (bigWigChromsFilename == None):
			bigWigChromsFilename = "{base}/temp/ucsc.hg19.chrom_lengths"
		if (bigWigUrl == None):
			usage("you have to give me a url for the bigwig file")

	trackId = "%s.clipped_breakpoints" % runName

	##########
	# perform filename substitution
	##########

	if   (basePath == None):       basePath = "."
	elif (basePath.endswith("/")): basePath = basePath[:-1]

	controlFilename = do_filename_substitutition(controlFilename)
	bamFilename     = do_filename_substitutition(bamFilename)
	chromsFilename  = do_filename_substitutition(chromsFilename)

	# track name

	trackName = do_filename_substitutition(trackName)

	trackFilename = trackName
	if (not trackFilename.endswith(".dat")): trackFilename += ".dat"

	# big wig name

	if (bigWigFilename != None):
		bigWigFilename = do_filename_substitutition(bigWigFilename)
		if ("{track}" in bigWigFilename):
			trackTemp = trackName
			if   (trackTemp.endswith(".gz")):  trackTemp = trackTemp[:-3]
			elif (trackTemp.endswith(".dat")): trackTemp = trackTemp[:-4]
			bigWigFilename = bigWigFilename.replace("{track}",trackTemp)
			if (bigWigFilename.endswith(".bw")): infoFilename = bigWigFilename[:-3] + ".info"
			else:                                infoFilename = bigWigFilename      + ".info"

	if (bigWigChromsFilename != None):
		bigWigChromsFilename = do_filename_substitutition(bigWigChromsFilename)

	if (bigWigUrl != None):
		bigWigTemp = bigWigFilename
		slashIx = bigWigTemp.rfind("/")
		if (slashIx >= 0): bigWigTemp = bigWigTemp[slashIx+1:]
		bigWigUrl = bigWigUrl.replace("{bigwig}",bigWigTemp)

	if (bigWigLink != None):
		bigWigSave = bigWigLink
		bigWigTemp = bigWigFilename
		slashIx = bigWigTemp.rfind("/")
		if (slashIx >= 0): bigWigTemp = bigWigTemp[slashIx+1:]
		bigWigLink = bigWigLink.replace("{bigwig}",bigWigTemp)

		infoTemp = infoFilename
		slashIx = infoTemp.rfind("/")
		if (slashIx >= 0): infoTemp = infoTemp[slashIx+1:]
		infoLink = bigWigSave.replace("{bigwig}",infoTemp)

	##########
	# get values from the control file
	##########

	clipThreshold = None
	windowSize    = None

	f = file(controlFilename,"rt")

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) >= 3), \
		       "not enough fields at control file line %d (%d, expected at least 3)" \
		     % (lineNumber,len(fields))
		assert (fields[1] == "="), \
		       "can't understand control file line %d:\n%s" \
		     % (lineNumber,line)

		(name,_,val) = fields[:3]
		if   (name == "clipThreshold"):                  clipThreshold = val
		elif (name == "clipped_breakpoints.windowSize"): windowSize    = val

	f.close()

	if (clipThreshold == None): assert(False), "control file lacks clipThreshold"
	if (windowSize    == None): assert(False), "control file lacks windowSize"

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

	if (dateInTrackname):
		print "today=`today {mmm}/{d}/{yyyy}`"

	# write commands describing the files the script will create

	print "echo \"will write track file to     %s\"" % trackFilename

	if (bigWigFilename != None):
		print "echo \"will write bigwig file to    %s\"" % bigWigFilename

	# write command(s) to create track file

	print
	print "echo \"=== creating track %s ===\"" % trackId

	commands =  []

	command  =  ["time samtools view %s" % bamFilename]
	commands += [command]

	command  =  ["awk '{ if (($6~/H/) || ($6~/S/)) print $0 }'"]
	commands += [command]

	command  =  ["filtered_sam_to_intervals"]
	if (isNameSorted):      command += ["--namesorted"]
	if (headLimit != None): command += ["--head=%s" % headLimit]
	command  += ["--prohibit:\"(CIGAR == *)\""]
	command  += ["--require:\" (RNEXT == =)\""]
	command  += ["--require:\" (UNCLIP > RLEN*%s)\"" % clipThreshold]
	command  += ["--report:\"  CLIPBRK\""]
	command  += ["--nonames"]
	command  += ["--progress=2M"]
	commands += [command]

	command  =  ["awk '{ print $1,$4,$4+1 }'"]
	commands += [command]

	command  =  ["genodsp"]
	command  += ["--chromosomes=%s" % chromsFilename]
	command  += ["--cliptochromosome --novalue --show:uncovered"]
	command  += ["= slidingsum W=%s" % windowSize]
	commands += [command]

	command  =  ["> %s" % trackFilename]
	commands += [command]

	print
	print commands_to_pipeline(commands)

	# write command(s) to convert track file to bigwig

	if (bigWigFilename != None):
		print
		print "echo \"=== converting track %s to bigwig ===\"" % trackId

		commands =  []
		command  =  ["time bedGraphToBigWig"]
		command  += [trackFilename]
		command  += [bigWigChromsFilename]
		command  += [bigWigFilename]
		commands += [command]

		print
		print commands_to_pipeline(commands)

		description = "local count of clipped breakpoints"
		if (dateInTrackname): description += " (${today})"

		commands =  []
		command  =  ["make_bigwig_info"]
		command  += ["--url=%s" % bigWigUrl]
		command  += ["--name=\"%s clipped breakpoints\"" % runName]
		command  += ["--desc=\"%s %s\"" % (runName,description)]
		command  += ["--autoscale=\"on\""]
		command  += ["--alwayszero=\"on\""]
		command  += ["--maxheight=\"128:30:20\""]
		command  += ["--color=120,20,230"]
		if (bigWigPosition != None): command  += ["--pos=\"%s\"" % bigWigPosition]
		command  += ["> %s" % infoFilename]
		commands += [command]

		print
		print commands_to_pipeline(commands)

		if (bigWigLink != None):
			print
			print "rm -f %s" % infoLink
			print "ln -s %s %s" % (infoFilename,infoLink)
			print "rm -f %s" % bigWigLink
			print "ln -s %s %s" % (bigWigFilename,bigWigLink)

			infoUrl = bigWigUrl
			slashIx = infoUrl.rfind("/")
			if (slashIx >= 0): infoUrl = infoUrl[:slashIx]
			infoTemp = infoFilename
			slashIx = infoTemp.rfind("/")
			if (slashIx >= 0): infoTemp = infoTemp[slashIx+1:]
			infoUrl = infoUrl + "/" + infoTemp
			print >>stderr, infoUrl
			print
			print "echo \"track URL is %s\"" % (infoUrl)


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

