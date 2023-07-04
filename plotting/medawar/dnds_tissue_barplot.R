##### dnds_tisue_barplot.R #########################################################################
# Plot summary statistics for beta_age vs dnds regression by medawarian trends
# 1) Barplot per tissue
# 2) Barplot per cell type

# Rscript ./dnds_barplot.R \
# -i /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/beta_dnds.csv \
# -t /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/regression/beta_dnds/dnds_tissue_regression.csv \
# -c /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/regression/beta_dnds/dnds_celltype_regression.csv


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
parser$add_argument('-b', '--beta', type = 'character', help = 'beta_dnds csv file');
parser$add_argument('-t', '--tissue', type = 'character', help = 'dnds tissue regression summary data csv file');
parser$add_argument('-c', '--celltype', type = 'character', help = 'dnds celltype regression summary data csv file');
args <- parser$parse_args();

#### READ DATA ####################################################################################
setwd('/global/scratch/users/pstein/Project-AgeExpressionConstraint/data/');
out.dir <- '/global/home/users/pstein/Project-AgeExpressionConstraint/figures/medawar/';

# Input file beta_age and dnds scores
beta.dnds <- read.csv(args$beta);
# Read in filtered version itself?
beta.dnds <- beta.dnds[beta.dnds$beta != 0,];
beta.dnds <- beta.dnds[!is.na(beta.dnds$dNdS),]; # 2531145

# Summary statistics of regression of dnds vs beta_age
df.tissue <- read.csv(args$tissue);
df.celltype <- read.csv(args$celltype);

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

# Genes per tissue
tissue <- sort(unique(beta.dnds$tissue));
n.gene.t <- data.frame(tissue);
n.gene.t$freq <- get.gene.count(beta.dnds, 'tissue');
df.tissue$n_genes <- paste0(n.gene.t$tissue, ' ', '(', as.character(n.gene.t$freq), ')');

# Genes per celltype
cell_type <- sort(unique(beta.dnds$cell_type));
n.gene.c <- data.frame(cell_type);
n.gene.c$freq <- get.gene.count(beta.dnds, 'cell_type');
df.celltype$n_genes <- paste0(n.gene.c$cell_type, ' ', '(', as.character(n.gene.c$freq), ')');

##### Barplot Tissue-specific Trends ##############################################################
plot.tissue <- ggplot(
  data = df.tissue, 
  aes(
    x = estimate, 
    y = reorder(n_genes, -estimate),
    fill = cond
  )
) +
  geom_bar(
    stat = "identity"
  ) + 
  labs(
    title = 'Tissue-specific Trends of Gene Expression in Mice', 
    x = expression(paste('Slope of dN/dS vs ',symbol('b')[Age])), 
    y = 'Tissue (n genes)'
  ) + 
  geom_errorbarh(
    aes(xmax = estimate + std.error, xmin = estimate - std.error, y = n_genes), 
    position=position_dodge(width = 0.75), 
    height = 0.5
  ) +
  geom_text(
    data = df.tissue,
    aes(x = estimate / abs(estimate) * (-0.02), y = n_genes, label = sig), 
    color = "black",
    nudge_y = -0.16,
    inherit.aes = FALSE,
    hjust= ifelse(df.tissue$estimate < 0, "inward", "outward"),
    size = 6
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
    filename = 'dnds_tissue_medawar_bar.pdf',
    plot = last_plot(),
    device = 'pdf',
    path = out.dir,
    width = 12,
    height = 8,
    units = 'in',
    dpi = 300,
    limitsize = TRUE
    );

##### Barplot Cell-specific Trends ################################################################
# Kinda a bad/ugly plot
plot.celltype <- ggplot(
  data = df.celltype, 
  aes(
    x = estimate, 
    y = reorder(n_genes, -estimate),
    fill = cond
    )
  ) +
  geom_bar(
    stat = "identity",
    color = 'white'
  ) + 
  labs(
    title = 'Cell-specific Trends of Gene Expression in Mice', 
    x = expression(paste('Slope of dN/dS vs ',symbol('b')[Age])), 
    y = 'Cell Type (n genes)'
  ) + 
  geom_errorbarh(
    aes(xmax = estimate + std.error, xmin = estimate - std.error, y = n_genes), 
    position=position_dodge(width = 0.75), 
    height = 0.5
  ) +
  geom_text(
    data = df.celltype,
    aes(x = estimate / abs(estimate) * (-0.001), y = n_genes, label = sig), 
    color = "black",
    nudge_y = -0.4,
    size = 3,
    inherit.aes = FALSE,
    hjust= ifelse(df.celltype$estimate < 0, "outward", "inward")
  ) + 
  scale_fill_manual(
    values = c('darkred', 'darkblue'),
    name = 'Trend'
  ) + 
  theme_bw() +
  theme(
    plot.title=element_text(hjust=0.5),
    axis.text = element_text(size=4.5),
    panel.background = element_blank()
  );
plot.celltype;

ggsave(
    filename = paste0('dnds_celltype_medawar_bar.pdf'),
    plot = last_plot(),
    device = 'pdf',
    path = out.dir,
    height = 8,
    width = 12,
    dpi = 300,
    limitsize = TRUE
    );
