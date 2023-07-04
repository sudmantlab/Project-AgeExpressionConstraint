import os
import json
import numpy as np

""" Generate a .json config file with the names of all the gene expression samples in the file folder. """

__author__ = "Philippa Steinberg"

cwd = os.getcwd()
output_dir = "../"

configfile = "/global/home/users/pstein/Project-AgeExpressionConstraint/analysis/expression/config.json"
with open(configfile, "r") as f:
        config = json.load(f)
        
def make_json():
    j_out = {} # initialize empty dict
    j_out["samples"] = {} # initialize empty sample (key) value (sample name) pair
    
    # directory with all mouse tisuse gene expression h5ad files
    sample_dir = config["TabulaMurisSenis"]
    for sample in os.listdir(sample_dir):
        sample_name = sample.split('-')[-1].split('.')[0]
        j_out["samples"][sample_name] = sample_dir + 'tabula-muris-senis-facs-processed-official-annotations-' + sample_name
    return j_out
FOUT = open(("{cwd}/" + "tissue_config.json").format(cwd = output_dir), "w")
FOUT.write(json.dumps(make_json(), indent = 4, separators = (",", ": ")))
FOUT.close()
