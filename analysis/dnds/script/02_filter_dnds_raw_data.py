import pandas as pd
import matplotlib as plt
import numpy as np
import os
from Bio import SeqIO

""" 
Process RELAX dNdS score output
1) Reformat data file 
2) Filter out genes with sequences < 100 bp (RELAX pipeline fix TBD)
"""

__author__ = "Philippa Steinberg"

configfile = "/global/home/users/pstein/Project-AgeExpressionConstraint/analysis/dnds/config.json"
with open(configfile, "r") as f:
        config = json.load(f)

# dNdS scores table (ensembl transcript ID version, dNdS score) from mouse-human-rat MSA
input_file = config["dnds_merged"]
# Directory with all gene MSA, aligned and formatted nucleotide codon MSA ending with "NT_CLEAN.fa"
sample_dir = config["NT_CLEAN.fa"]
# dNdS scores table (ensembl transcript ID version, ensembl transcript ID, dNdS score)
output_file = 'ensmust_dnds_filtered_100bp.csv'

# Read in data and rename indices
dnds_scores = pd.read_csv(input_file, sep ='\t', header=None, names=['ENSMUST.v', 'dNdS'])
dnds_scores['ENSMUST'] = list(dnds_scores['ENSMUST.v'].str.split('.').str[0])

# Remove rows with dNdS = 0 (no substitutions)
dnds_scores = dnds_scores.set_index(['ENSMUST', 'ENSMUST.v'])
dnds_scores= dnds_scores.round(6).loc[~(dnds_scores==0).all(axis=1)]

# Function to create list of MSA with gene seq length at least 100 bp
def filter_100_bp(directory):

    ensmust_list = []
    for sample in os.listdir(directory):
        # Check for aligned and formatted MSA samples
        if sample.endswith("NT_CLEAN.fa"):
            filepath = os.path.join(directory, sample)
            records = SeqIO.parse(filepath, "fasta")
            # check only first seq in MSA (all taxa should have same len)
            first_record = next(records)
            # check that seq > 100bp
            if len(first_record.seq) > 100:
                # print(f"{sample}'s DNA sequence is longer than 100bp")
                # get only exon ID
                exon_name = sample.split("_")[0]
                ensmust_list.append(exon_name)
    return ensmust_list
    
ensmust_100 = filter_100_bp(sample_dir)
dnds_scores = dnds_scores.reset_index()
dnds_scores_100 = dnds_scores.loc[dnds_scores['ENSMUST.v'].isin(ensmust_100)]
dnds_scores_100 = dnds_scores_100.set_index(['ENSMUST', 'ENSMUST.v'])

stats = dnds_scores_100.describe()
print(stats)
dnds_scores_100.to_csv(output_file)
