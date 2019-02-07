#!/usr/bin/python
"""
Report today's date.
"""

from sys  import argv,stdin,stderr,exit
from time import strftime


def usage(s=None):
	message = """
usage: today <format>
  <format> is a format string suitable for python's strftime function.
  --head=<number>        limit the number of input lines
  --progress=<number>    periodically report how many lines we've read"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	if (len(argv) == 1): formats = ["{mmm}/{d}/{yyyy}"]
	else:                formats = argv[1:]

	for format in formats:
		format = format.replace("{yyyy}",    "%Y")	# year

		format = format.replace("{month}",   "%B")	# month
		format = format.replace("{mmm}",     "%b")	# month as three letters
		format = format.replace("{mm}",      "%m")	# month as two digits
		format = format.replace("{m}",       "~%m")	# month as number

		format = format.replace("{weekday}", "%A")	# weekday
		format = format.replace("{www}",     "%a")	# weekday as three letters

		format = format.replace("{day}",     "~%d")	# day of month
		format = format.replace("{dd}",      "%d")	# day of month as two digits
		format = format.replace("{d}",       "~%d")	# day of month

		format = format.replace("{hms}",     "%I:%M:%S %p")
		format = format.replace("{hms24}",   "%H:%M:%S")

		format = format.replace("{hour}",    "%I")
		format = format.replace("{hr}",      "%I")
		format = format.replace("{hour12}",  "%I")
		format = format.replace("{hr12}",    "%I")
		format = format.replace("{hour24}",  "%H")
		format = format.replace("{hr24}",    "%H")

		format = format.replace("{minute}",  "%M")
		format = format.replace("{min}",     "%M")

		format = format.replace("{second}",  "%S")
		format = format.replace("{sec}",     "%S")

		format = format.replace("{am}",      "%p")
		format = format.replace("{pm}",      "%p")

		today = strftime(format)
		today = today.replace("~0","")
		today = today.replace("~","")
		print today


if __name__ == "__main__": main()
