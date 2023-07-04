## Guide
Calculate normalized beta_age values per gene for each tissue.
- config.json contains paths to directories and cell_types which were excluded from the tissue analysis 
- Run beta_age.ipynb to calculate beta_age per tissue (modify counter and run interactively)
- Run generate_config_file.json to generate tissue_config.json (useful if running beta_age.ipnyb with snakemake)
- tissue_config.json contains path to 23 TabulaMurisSenis (facs) gene expression tissues as h5ad files
