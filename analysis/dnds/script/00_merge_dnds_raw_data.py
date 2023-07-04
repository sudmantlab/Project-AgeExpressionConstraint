import json
import os

""" concatenate individual ensmust_dnds.txt RELAX output files to one txt file """

__author__ = "Philippa Steinberg"

configfile = "/global/home/users/pstein/Project-AgeExpressionConstraint/analysis/dnds/config.json"
with open(configfile, "r") as f:
        config = json.load(f)

dnds_dir = config["dnds.txt"]
output_file = 'ensmust_dnds_merged.txt'

with open(output_file, 'w') as f_out:
	for filename in os.listdir(dnds_dir):
		if filename.endswith(".txt"):
			ensmust = filename.split("_")[0]
			filepath = os.path.join(dnds_dir, filename)
			with open(filepath,'r') as f_in:
				file_contents = f_in.readline()
				f_out.write(ensmust + '\t')
				f_out.write(file_contents + '\n')