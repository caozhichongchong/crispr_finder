# crispr_finder
Code for paper "CRISPR spacer acquisition is a rare event in human gut microbiome"
## Code for data summarization and visualization
* mathematical model: crispr_model.ipynb
* summarizing spacer rates in metagenomic data: crispr_MG_rate.ipynb
* summarizing spacer rates in WGS data: crispr_WGS_rate.ipynb
* summarizing HGT of spacers in B.longum species: crispr_HGT.ipynb
## Code for sequencing data analysis - pipeline in crispr.py
* running cirsprcasefinder for WGS data: runCRISPRCasFinder.py
* merging and summarizzing crispr spacers + direct repeats: merge_spacer_filter.py
* detecting spacer acquisition in WGS data using maximum likelihood model: spacer_parsimony.py and spacer_parsimony_sum.py
* annotating spacers: annotate_spacer.py
* constructing the structure of CRISPR systems: spacer_structure.py
* extracting spacers + direct repeats from a specific species: DR_extract.py
* mapping direct repeats in metagenomes to WGS data: DR_MGmapping.py
## Data
* containing data required to run jupyter notebooks
