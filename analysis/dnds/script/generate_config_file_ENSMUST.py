import os
import json
import numpy as np

"""
Generate a .json config file with the names of all the samples in the ENSMUST_fa file folder that correspond to the first (cannonical) exon. 
"""

__author__ = "Philippa Steinberg"

configfile = "/global/home/users/pstein/Project-AgeExpressionConstraint/analysis/dnds/config.json"
with open(configfile, "r") as f:
        config = json.load(f)

cwd = os.getcwd()
output_dir = "/global/scratch/users/pstein/Project-AgeExpressionConstraint/data/dnds/"

def make_json():
	j_out = {} # initialize empty dict
	j_out["samples"] = {} # initialize empty sample (key) value (sample name) pair
	
	# directory with all fasta files by ENMUST and exon for 5 taxa of interest
	sample_dir = config["ENSMUST_fa"]
	for sample in os.listdir(sample_dir):
		# select only files with first (canonical) exon
		exon = sample.split("_")[1]
		if exon == "1":
			# remove "_ID.fa" extension
			sample_name = sample[:-6]
			j_out["samples"][sample_name] = sample_dir + sample_name
	return j_out

FOUT = open(("{cwd}/" + "ENSMUST_all.json").format(cwd = output_dir), "w")
FOUT.write(json.dumps(make_json(), indent = 4, separators = (",", ": ")))
FOUT.close()
