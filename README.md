# Age vs Expression vs Constraint 
## Project Overview

This project is a continuation of project **gene_expression_aging** https://github.com/sudmantlab/gene_expression_aging]. 
Contrasting selective constraint (phyloP) and selection (dNdS) with age-associated gene expression for mice at the tissue and celltype level. 
Includes dNdS calculation of mouse genome, Medawarian trend and Hallmark (GSEA) analysis.

## Background
The foundation of my research question rests on Medawar’s hypothesis. Under this theory, we expect reduced purifying selection and thus increased accumulation of deleterious mutations for genes expressed late in life. In contrast, functionally important genes expressed early in life are expected to be more conserved (Charlesworth, 2000). The paper by Charlesworth investigated the rate of fixation of deleterious mutations across age and confirms Medawar’s hypothesis for humans and mice on the organism level. The paper the Sudmant lab (Yamamoto and Chung et al. 2022) contrasts tissue-specific gene expression throughout age with selective constraint measured as being loss-of-function intolerant (pLI). Further studies analyze aging on a cell-specific level in mice to distinguish if changes in gene expression are due to a change in cell quantity over time or cell-specific gene expression with age (Schaum et al. 2020). This project is an analysis of the tissue- and cell-specific trends of age-associated gene expression. While Yamamoto and Chung use pLI, this project uses PhyloP and dNdS scores as metrics of constraint. 


## Folder Structure

```
analysis/
|-- README.md
|-- expression/
	|-- ..
|-- dnds/
	|-- ..
|-- phyloP/
	|-- ..
|-- age_expression_dnds_merge.R
|-- age_expression_phyloP_merge.R
|-- beta_dnds_regression.R
|-- beta_phyloP_regression.R
```
- **expression**: Code for calculating `beta_age` per gene for TabulaMurisSenis gene expression data.
- **dnds**: Code for calculating `dNdS` scores per mouse gene using RELAX in a Snakemake pipeline.
- **phyloP**: Code for calculating average `phyloP` per exon for mouse phyloP data.

```
plotting/
|-- README.md
|-- correlation/
	|-- dnds_phyloP.R
	|-- dnds_scatterplot.R
	|-- phyloP_scatterplot.R
|-- medawar/
	|-- dnds_tissue_barplot.R
	|-- dnds_celltype_barplot.R
	|-- phyloP_tissue_barplot.R
	|-- dnds_celltype_barplot.R
|-- gsea/
	|-- dnds_gsea.R
	|-- phyloP_gsea.R

```
- **correlation**: Code for correlating metrics of selection and for correlating `beta_age` with selection.
- **medawar**: Code for creating per tissue and per celltype medawarian trends of `beta_age` vs selection.
- **gsea**: Code for hallmark analysis by young and aged.


## Data Structure

All project input and output files are on savio scratch directory: `/global/scratch/users/pstein/Project-AgeExpressionConstraint/`


```
data-raw/
|-- expression/
	|-- TabulaMurisSenis_h5ad/
		|-- tabula-muris-senis-facs-processed-official-annotations-{tissue}.h5ad
|-- phyloP/
	|-- mm10.60way.phyloP60way.bw
	|-- knownGene.txt
|-- dnds/
	|-- dnds_paper_vidal.xls
	|-- knownGene.exonNuc.fa.gz
	|-- species_TimeTree_renamed_3.nwk
	|-- taxa_names.lst
	|-- intermediate/ (intermediate.zip)
		|-- {ensmust.v}_TAXA.fa
		|-- {ensmust.v}_NT.fa
		|-- {ensmust.v}_AA.fa
		|-- {ensmust.v}_NT_CLEAN.fa
		|-- {ensmust.v}_AA_CLEAN.fa
	|-- relax/ (relax.zip)
		|-- {ensmust.v}_NT_RELAX.json
|-- id_conversion/
	|-- mgi_to_ensembl_transcript_id.csv
```
- **TabulaMurisSenis_h5ad**: https://figshare.com/articles/dataset/Tabula_Muris_Senis_Data_Objects/12654728?file=23872460
- **phyloP**: http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phyloP60way/, https://hgdownload.soe.ucsc.edu/gbdb/mm10/
- **dnds/dnds_paper.csv**: https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-599#MOESM1
- **dnds/knownGene.exonNuc.fa**: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/multiz60way/alignments/ and separated by gene, subset by taxa using `split_ensmust.txt`, aligned, exported and calculated dNdS scores with relax.


```
data/
|-- beta_dnds.csv
|-- beta_phyloP_mean.csv
|-- beta_phyloP_first.csv
|-- beta_age/
	|-- {tissue}_age_expression_normed.csv
	|-- beta_age_normalized.csv
|-- phyloP/
	|-- mus_exon_phyloP_mean.csv
  	|-- mus_exon_phyloP_first.csv
|-- dnds/
	|-- ENSMUST_all.json
	|-- ensmust_dnds_merged.txt
	|-- ensmust_dnds_processed.csv
	|-- ensmust_dnds_filtered_100bp.csv
	|-- data/ (data.zip)
|-- regression/
	|-- beta_dnds/
		|-- dnds_tissue_regression.csv
		|-- dnds_celltype_regression.csv
		|-- dnds_tissue_celltype_regression.csv
	|-- beta_phyloP/
		|-- phyloP_tissue_regression.csv
		|-- phyloP_celltype_regression.csv
		|-- phyloP_tissue_celltype_regression.csv
	
```
- **beta_dnds** : All `beta_age` and `dnds` (len > 100bp) values merged by gene. 
- **phyloP_dnds** : All `beta_age` and `phyloP` scores merged by gene (version of mean phyloP and only first exon phyloP score). 

## Data Processing Comments
- Age categories 3m (young), 24m (aged)
- Gene expression values were z-score normalized
- PhyloP scores were calculated as mean across all exons and only the first exon. 
- Mouse dNdS scores from paper seem low and thus not used.
- dNdS scores calculated with RELAX (90% complete). Some short genes have dNdS >>> scores. View issue posted on HyPhy developers repo and their reply indicating the fix TBA https://github.com/veg/hyphy/issues/1614

## Figures 
Figures generated by analysis are in folder `/global/home/users/pstein/Project-AgeExpressionConstraint/figures/`.

Figure 1. `./correlation`
- a. correlate phyloP vs dnds (paper)
- b. correlate phylop vs dnds (relax)
- c. correlate dnds (paper) vs dnds (relax)

Figure 2. `./correlation`
- a. scatterplot celltype specific trends example tissue dnds
- b. barplot celltype specific trends example tissue dnds
- c. scatterplot celltype specific trends example tissue phyloP
- d. barplot celltype specific trends example tissue phyloP

Figure 3. `./medawar/`
- a. barplot tissue specific trends dnds
- b. barplot tissue specific trends phyloP

Figure 4. `./gsea`
- a. gsea dnds
- b. gsea phyloP


## Author information

- savio: pstein
- uc berkeley email: p.stein@berkeley.edu (active)
- uw email: pstein@uw.edu (active)

Template credit for `generate_config_file.py` to Alma Hallgren and template credit for [gsea.R] to Peter Sudmant.