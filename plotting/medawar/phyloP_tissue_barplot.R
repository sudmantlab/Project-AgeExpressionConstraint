##### phyloP_tissue_barplot.R #######################################################################
# Plot summary statistics for beta_age vs phyloP regression by medawarian trends
# 1) Barplot per tissue
# 2) Barplot per cell type

# Rscript ./phyloP_tissue_barplot.R \
# -i /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/beta_phyloP_mean.csv \
# -t /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/regression/beta_phyloP/phyloP_mean_tissue_regression.csv \
# -c /global/scratch/users/pstein/Project-AgeExpressionConstraint/data/regression/beta_phyloP/phyloP_mean_celltype_regression.csv

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
parser$add_argument('-b', '--beta', type = 'character', help = 'beta_phyloP (first or mean) csv file');
parser$add_argument('-t', '--tissue', type = 'character', help = 'phyloP tissue regression summary data csv file');
parser$add_argument('-c', '--celltype', type = 'character', help = 'phyloP celltype regression summary data csv file');
args <- parser$parse_args();

#### READ DATA ####################################################################################
setwd('/global/scratch/users/pstein/Project-AgeExpressionConstraint/data/');
out.dir <- '/global/home/users/pstein/Project-AgeExpressionConstraint/figures/medawar/';

# Input file beta_age and phyloP scores
beta.phyloP <- read.csv(args$beta);

# Read in filtered version itself?
beta.phyloP <- beta.phyloP[beta.phyloP$beta != 0,];
# If beta_phyloP_first.csv
# beta.phyloP <- beta.phyloP[!is.na(beta.phyloP$Mouse.first.exon.phyloP),];
# If beta_phyloP_mean.csv
beta.phyloP <- beta.phyloP[!is.na(beta.phyloP$Mouse.exon.phyloP.mean),]; # 13057492

df.tissue <- read.csv(args$tissue);
df.celltype <- read.csv(args$celltype);

##### Add Num Genes to DFs ########################################################################
# function to get gene count
get.gene.count <- function(df, group) {
    gene_count <- list()
    for (t in sort(unique(df[[group]]))) {
        data <- df[df[[group]] == t,]
        n_gene <- length(unique(data$mgi_symbol))
        gene_count <- append(gene_count, n_gene)
        }
    return(gene_count)
    };

# Genes per tissue
tissue <- sort(unique(beta.phyloP$tissue));
n.gene.t <- data.frame(tissue);
n.gene.t$freq <- get.gene.count(beta.phyloP, 'tissue');
df.tissue$n_genes <- paste0(n.gene.t$tissue, ' ', '(', as.character(n.gene.t$freq), ')');

# Genes per celltype
cell_type <- sort(unique(beta.phyloP$cell_type));
n.gene.c <- data.frame(cell_type);
n.gene.c$freq <- get.gene.count(beta.phyloP, 'cell_type');
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
    x = expression(paste('Slope of phyloP vs ',symbol('b')[Age])), 
    y = 'Tissue (n genes)'
  ) + 
  geom_errorbarh(
    aes(xmax = estimate + std.error, xmin = estimate - std.error, y = n_genes), 
    position=position_dodge(width = 0.75), 
    height = 0.5
  ) +
  geom_text(
    data = df.tissue,
    aes(x = estimate / abs(estimate) * (-0.05), y = n_genes, label = sig), 
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
    filename = 'phyloP_tissue_medawar_bar.pdf',
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
    aes(x = estimate / abs(estimate) * (-0.005), y = n_genes, label = sig), 
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
    filename = 'phyloP_celltype_medawar_bar.pdf',
    plot = last_plot(),
    device = 'pdf',
    path = out.dir,
    height = 8,
    width = 12,
    dpi = 300,
    limitsize = TRUE
    );
