#### age_expression_dnds_merge.R ##################################################################
# Merge beta_age (expression) values with dNdS scores via the MGI to ENSMUST translation table
# Rscript ./age_expression_dnds_merge.R \
# -d /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/dnds/ensmust_dnds_filtered_100bp.csv \
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
parser$add_argument('-d', '--dNdS', type = 'character', help = 'dNdS scores csv file (ENSMUST, dNdS)');
parser$add_argument('-b', '--beta', type = 'character', help = 'beta_age (expression) scores csv file');
parser$add_argument('-t', '--translation', type = 'character', help = 'mgi symbol to ensembl transcript ID translation csv file');
args <- parser$parse_args();

#### READ DATA ####################################################################################
setwd('/global/scratch/users/pstein/Project-AgeExpressionConstraint/data/');
# Mouse dNdS data
dnds.scores <- read.csv(args$dNdS);

# Beta age values for all tissues per gene
age.expression.scores <- read.csv(args$beta);
# age.expression.scores <- read.csv('./beta_age/beta_age_normalized.csv');

# MGI symbol to ENSMUST mapping table
translation.ids <- read.csv(args$translation);
# translation.ids <- read.csv('./id_conversion/mgi_to_ensembl_transcript_id.csv');

##### MERGE #######################################################################################
# Map mgi symbol to translation
mapping <- merge(
    age.expression.scores, 
    translation.ids, 
    by.x = 'gene',
    by.y = 'mgi_symbol'
    );

# Map ensmust with dnds 
expression.dnds <- merge(
    mapping,
    dnds.scores,
    by.x = 'ensembl_transcript_id',
    by.y = 'ENSMUST'
    );

##### FORMAT #################################################################################
expression.dnds.df <- expression.dnds[, c(
    'gene',
    'ensembl_gene_id',
    'ensembl_transcript_id',
    'ENSMUST.v',
    'tissue',
    'cell_type',
    'rho',
    'beta', 
    'p_value',
    'dNdS'
    )];

# Save output
write.csv(expression.dnds.df, "./beta_dnds.csv", row.names = FALSE);
