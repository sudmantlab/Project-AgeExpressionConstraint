##### dnds_celltype_barplot.R #####################################################################
# Plot summary statistics for beta_age vs dnds regression by medawarian trends
# 1) Barplot per cell type for each tissue

# Rscript ./dnds_celltype_barplot.R \
# -b /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/beta_dnds.csv \
# -s /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/regression/beta_dnds/dnds_celltype_tissue_regression.csv \
# -o /global/home/users/pstein/Project-AgeExpressionConstraint/figures/medawar/ \
# -t Liver

##### Preamble ####################################################################################
# Libraries
library(dplyr);
library(plyr);
library(tidyr);
library(tidyverse);
library(argparse);
library(rjson);
library(R.utils);

#### OBTAIN COMMAND LINE ARGUMENTS ################################################################
parser <- ArgumentParser();
parser$add_argument('-b', '--beta', type = 'character', help = 'beta_dnds csv file');
parser$add_argument('s', '--summary', type = 'character', help = 'regression summary file csv');
parser$add_argument('-o', '--outdir', type = 'character', help = 'output directory of figure');
parser$add_argument('t', '--tissue', type = 'character', help = 'tissue from config.json');
args <- parser$parse_args();

#### READ DATA ####################################################################################
setwd('/global/scratch/users/pstein/Project-AgeExpressionConstraint/data/');
out.dir <- '/global/home/users/pstein/Project-AgeExpressionConstraint/figures/medawar/';

# Input file beta_age and dnds scores
beta.dnds <- read.csv(args$input);

# Read in filtered version itself?
beta.dnds <- beta.dnds[beta.dnds$beta != 0,];
beta.dnds <- beta.dnds[!is.na(beta.dnds$dNdS),]; # 2531145

# Summary statistics of regression of dnds vs beta_age
df.celltype.tissue <- read.csv(args$celltypetissue);

##### Add Num Genes to DFs ########################################################################
# function to get gene count
get.gene.count <- function(df, group) {
    gene_count <- list()
    for (t in sort(unique(df[[group]]))){
        data <- df[df[[group]] == t,]
        n_gene <- length(unique(data$gene))
        gene_count <- append(gene_count, n_gene)
        }
    return(gene_count)
    };

##### Plot Cell Types for Specific Tissue #########################################################
# Subset to tissue of interest
tissue <- args$tissue;
tissue_data <- df.celltype.tissue[df.celltype.tissue$tissue == tissue,];
tissue_data <-tissue_data[order(tissue_data$cell_type),];
beta_data <- beta.dnds[beta.dnds$tissue == tissue,]; # subset to only tissue of interest

# Genes per celltype for substetted tissue
cell_type <- sort(unique(tissue_data$cell_type));
n.gene.c <- data.frame(cell_type);
n.gene.c$freq <- get.gene.count(beta_data, 'cell_type');
tissue_data$n_genes <- paste0(n.gene.c$cell_type, ' ', '(', as.character(n.gene.c$freq), ')');

### Barplot celltype-specific trends for one tissue ###############################################
# Barplot (Add tissue name manually in title, modify spacing of "***" manually)
plot.tissue <- ggplot(
  data = tissue_data, 
  aes(
    x = estimate, 
    y = reorder(n_genes, -estimate), 
    fill = cond)
  ) +
  geom_bar(
    stat = "identity"
  ) + 
  labs(
    title = paste0('Cell-specific Trends of Gene Expression in ', tissue), 
    x = expression(paste('Slope of dN/dS vs ',symbol('b')[Age])), 
    y = 'Cell Type (n genes)'
  ) + 
  geom_errorbarh(
    aes(xmax = estimate + std.error, xmin = estimate - std.error, y = n_genes), 
    position=position_dodge(width = 0.75), 
    height = 0.5
  ) + 
  geom_text(
    data = tissue_data,
    aes(x = estimate / abs(estimate) * (-0.03), y = n_genes, label = sig), 
    color = "black", size = 6,
    nudge_y = -0.02,
    inherit.aes = FALSE,
    hjust= ifelse(tissue_data$estimate < 0, "inward", "outward")
  ) + 
  scale_fill_manual(
    values = c('darkred', 'darkblue'),
    name = "Trend"
  ) + 
  theme_bw() +
  theme(
    plot.title=element_text(hjust=0.5)
  );
plot.tissue;

ggsave(
    filename = paste0(tissue, '_dnds_medawar_bar.pdf'),
    plot = last_plot(),
    device = 'pdf',
    path = out.dir,
    height = 8,
    width = 12,
    dpi = 300,
    limitsize = TRUE
    );
