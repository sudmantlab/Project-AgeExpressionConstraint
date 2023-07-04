##### beta_dnds_regression.R ######################################################################
# Calculate summary statistics for beta_age vs dnds regression
# 1) Table per tissue
# 2) Table per cell type
# 3) Table per cell type for each tissue

# Rscript ./beta_dnds_regression.R \
# -i /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/beta_dnds.csv \
# -o /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/regression/beta_dnds

##### Preamble ####################################################################################
# Libraries
library(broom);
library(dplyr);
library(plyr);
library(tidyr);
library(tidyverse);
library(argparse);
library(rjson);
library(R.utils);

#### OBTAIN COMMAND LINE ARGUMENTS ################################################################
parser <- ArgumentParser();
parser$add_argument('-i', '--input', type = 'character', help = 'beta_dnds csv file');
parser$add_argument('-o', '--output', type = 'character', help = 'regression output directory');
args <- parser$parse_args();

#### READ DATA ####################################################################################
setwd('/global/scratch/users/pstein/Project-AgeExpressionConstraint/data/');

# Input file beta_age and dnds scores
#beta.dnds <- read.csv('./beta_dnds.csv');
beta.dnds <- read.csv(args$input);

# Output directort 
#out.dir <- '/global/scratch/users/pstein/Project-AgeExpressionConstraint/data/regression/beta_dnds/';
out.dir <- args$output;

##### Tissue Stats DF #############################################################################
# Create and save a df with regression (dnds vs beta_age) statistics per tissue
df_tissue <- data.frame();
for (tissue in unique(beta.dnds$tissue)){
    print(tissue)
    data <- beta.dnds[beta.dnds$tissue == tissue,];
    stats <- data %>% do(
        tidy(lm(data$dNdS ~ data$beta))
        )
    stats$tissue <- tissue
    print(stats)
    df <- data.frame(stats)
    df_tissue <- rbind(df_tissue,df[2,])
    };

##### Cell Type Stats DF ##########################################################################
# Create and save a df with regression (dnds vs beta_age) statistics per cell type
df_celltype <- data.frame();
for (cell in unique(beta.dnds$cell_type)){
    print(cell)
    data <- beta.dnds[beta.dnds$cell_type == cell,];
    stats <- data %>% do(
        tidy(lm(data$dNdS ~ data$beta))
        )
    stats$cell_type <- cell
    print(stats)
    df <- data.frame(stats)
    df_celltype <- rbind(df_celltype, df[2,])
    };

##### Cell Type per Tissue Stats DF ###############################################################
# Create and save a df with regression (dnds vs beta_age) statistics per tissue and cell type
df_celltype_tissue <- data.frame();
for (tissue in unique(beta.dnds$tissue)){
    print(tissue)
    tissue_data <- beta.dnds[beta.dnds$tissue == tissue,];
    for (cell in unique(tissue_data$cell_type)){
        print(cell)
        cell_data <- tissue_data[tissue_data$cell_type == cell,];
        stats <- cell_data %>% do(
            tidy(lm(cell_data$dNdS ~ cell_data$beta))
            )
        stats$tissue <- tissue
        stats$cell_type <- cell
        df <- data.frame(stats)
        df_celltype_tissue <- rbind(df_celltype_tissue, df[2,])
        };
    };

##### Add Medawarian Trends and P-values ##########################################################
# Could make a common fn so I do not need to repeat it each time
# Add medawarian trends to dfs
df_tissue <- df_tissue %>%
    mutate(cond = ifelse(
        estimate > 0, 
        'non-medawarian',
        'medawarian' 
        )
    );

df_celltype <- df_celltype %>%
    mutate(cond = ifelse(
        estimate > 0, 
        'non-medawarian',
        'medawarian' 
        )
    );

df_celltype_tissue <- df_celltype_tissue %>%
    mutate(cond = ifelse(
        estimate > 0, 
        'non-medawarian',
        'medawarian' 
        )
    );

# Add p-value significance to DFs
df_tissue <- df_tissue %>% mutate(
    sig = ifelse(p.value > 0.05, '',
                 ifelse(p.value > 0.01, '*',
                        ifelse(p.value > 0.001, '**', '***')
                 )
    )
);

df_celltype <- df_celltype %>% mutate(
    sig = ifelse(p.value > 0.05, '',
               ifelse(p.value > 0.01, '*',
                      ifelse(p.value > 0.001, '**', '***')
               )
    )
);

df_celltype_tissue <- df_celltype_tissue %>% mutate(
    sig = ifelse(p.value > 0.05, '',
               ifelse(p.value > 0.01, '*',
                      ifelse(p.value > 0.001, '**', '***')
               )
    )
);

##### SAVE AS DFS  ################################################################################
save <- write.csv(df_tissue, paste0(out.dir, 'dnds_tissue_regression.csv'), row.names = FALSE);
save <- write.csv(df_celltype, paste0(out.dir, 'dnds_celltype_regression.csv'), row.names = FALSE);
save <- write.csv(df_celltype_tissue, paste0(out.dir, 'dnds_celltype_tissue_regression.csv'), row.names = FALSE);
