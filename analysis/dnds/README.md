## Guide

Need to adjust paths according to config.json.
Need to unzip processed files in data-raw/dnds/ for further use. 

1. `split_ensmust.txt` contains code for splitting whole genome MSA into per gene MSA fasta files in `ENSMUST_fa` folder.
2. Run `generate_config_file_ENSMUST.py` to Generate ENSMUST_all.json config with all gene fa files from `ENSMUST_fa` folder. 
3. Run `Snakefile_filter` to subset MSA for taxa of interest (mouse, rat, human)
4. Species tree for dNdS calculation in tree_files
5. Run `Snakefile_full` to align, export, calculate dNdS (RELAX) and extract dNdS scores per gene (cannonical exon)
6. Run `00_merge_dnds_raw_data.py` to generate dnds output table.
7. Run `01_process_dnds_raw_data.py` to format dnds output table.
8. Run `02_filter_dnds_raw_data.py` to filter for genes > 100 bp. 
