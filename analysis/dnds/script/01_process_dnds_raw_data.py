import pandas as pd
import matplotlib as plt
import numpy as np
import os

"""
Process RELAX dNdS score output
1) Reformat data file 
2) Filter out all dNdS = 0 (no substitutions)
"""

__author__ = "Philippa Steinberg"

configfile = "/global/home/users/pstein/Project-AgeExpressionConstraint/analysis/dnds/config.json"
with open(configfile, "r") as f:
        config = json.load(f)

# dNdS scores table (ensembl transcript ID version, dNdS score) from mouse-human-rat MSA
input_file = config["dnds_merged"]
# dNdS scores table (ensembl transcript ID version, ensembl transcript ID, dNdS score)
output_file = 'ensmust_dnds_processed.csv'

# Read in data and rename indices
dnds_scores = pd.read_csv(input_file, sep ='\t', header=None, names=['ENSMUST.v', 'dNdS'])
dnds_scores['ENSMUST'] = list(dnds_scores['ENSMUST.v'].str.split('.').str[0])
dnds_scores = dnds_scores.set_index(['ENSMUST', 'ENSMUST.v'])

# Remove rows with dNdS = 0 (no substitutions)
dnds_scores = dnds_scores.round(6).loc[~(dnds_scores==0).all(axis=1)]

stats = dnds_scores.describe()
print(stats)
dnds_scores.to_csv(output_file)


