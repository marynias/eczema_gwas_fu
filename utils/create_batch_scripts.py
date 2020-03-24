#! /usr/bin/env python
import os, sys, re

def main(script_template, script_output_name, **kwargs):
	out = open(script_output_name, 'w')
	with open(script_template, 'r') as st_h:
		for line in st_h:
			if line.startswith("%%"):
				lines = line.strip().split()
				arguments = lines[1].split("=")
				if arguments[0] in kwargs:
					if arguments[0] == "optional":
						out.write(arguments[0] + "=" + '"' + kwargs[arguments[0]] + '"' + "\n")
					else:
						out.write(arguments[0] + "=" + kwargs[arguments[0]] + "\n")
				#else:
				#	raise ValueError ("Missing value for argument %s" % arguments[0])
			else:
				out.write(line)
	out.close()        

if __name__=='__main__':
    main(sys.argv[1], # foo
         sys.argv[2], # bar
         **dict(arg.split('=') for arg in sys.argv[3:])) # kwargs

