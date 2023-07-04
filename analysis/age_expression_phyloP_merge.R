##### age_expression_phyloP_merge.R ###############################################################
# Merge beta_age (expression) values with phyloP scores via the MGI to ENSMUST translation table
# Rscript ./age_expression_phyloP_merge.R \
# -d /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/phyloP/mus_exon_phyloP_mean.csv \
# -b /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/beta_age/beta_age_normalized.csv \
# -t /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/id_conversion/mgi_to_ensembl_transcript_id.csv

##### Preamble ####################################################################################
# install.packages('argparse');
# install.packages('rjson');
# install.packages('R.utils');

library(argparse);
library(rjson);
library(R.utils);

#### OBTAIN COMMAND LINE ARGUMENTS ################################################################
parser <- ArgumentParser();
parser$add_argument('-d', '--phyloP', type = 'character', help = 'phyloP scores csv file (Transcript.stable.ID.version, Mouse.exon.phyloP.mean/Mouse.exon.phyloP.first)');
parser$add_argument('-b', '--beta', type = 'character', help = 'beta_age (expression) scores csv file');
parser$add_argument('-t', '--translation', type = 'character', help = 'mgi symbol to ensembl transcript ID translation csv file');
args <- parser$parse_args();

#### READ DATA ####################################################################################
setwd('/global/scratch/users/pstein/Project-AgeExpressionConstraint/data/');

# Mouse PhyloP data
phyloP.scores <- read.csv(args$phyloP);
# phyloP.scores <- read.csv('./phyloP/mus_exon_phyloP_first.csv');

# Beta age values for all tissues per gene
age.expression.scores <- read.csv(args$beta);
# age.expression.scores <- read.csv('./beta_age/beta_age_normalized.csv');

# MGI symbol to ENSMUST mapping table
translation.ids <- read.csv(args$translation);
# translation.ids <- read.csv('./id_conversion/mgi_to_ensembl_transcript_id.csv');

##### Merge data frames ###########################################################################
# Merge phyloP scores to ENSMUST in mapping table
mapping <- merge(
    phyloP.scores, 
    translation.ids, 
    by.x = 'Transcript.stable.ID', 
    by.y = 'ensembl_transcript_id'
    );

# Merge MGI symbol with phyloP data
expression.phyloP <- merge(
  mapping, 
    age.expression.scores, 
    by.x = 'mgi_symbol', 
    by.y = 'gene'
    );

# Save table of merged beta and phyloP scores by gene
save <- write.csv(expression.phyloP, "./beta_phyloP.csv", row.names = FALSE);
