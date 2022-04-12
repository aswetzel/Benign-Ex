# Benign-Ex
Benign-Ex is designed to classify regions of the genome as “benign” from the perspective of copy number status.  Benign-Ex relies on two assumptions.  First, that CNVs commonly classified clinically as benign or likely benign occur in regions of the genome that are benign from a copy number status perspective.  Second, that regions of the genome that exhibit a high prevalence of copy number variation in the general population are benign from a copy number perspective. 

Benign-Ex can be run in two ways:

1)	Based on a single set of parameters provided by the user, Benign-Ex will identify normal copy number variable (“benign”) regions of the human genome from frequency-based and/or classification-based datasets. 

2)	Based on multiple sets of parameters provided by the user, Benign-Ex will identify the optimal set of parameters for identifying normal copy number variable (“benign”) regions of the human genome from frequency-based and/or classification-based datasets. The optimal set of parameters is determined by comparing the amount of overlap between the Benign-Ex identified “benign” regions with a set or sets of genomic regions known (or highly likely) to be pathogenic from a copy number perspective. 

# Requirements
Benign-Ex requires Python (2.7.X) and R (>=3.6.0) in a Linux operating system to run. Python3 is currently not supported. A list of the required Python modules and R libraries is listed below. Benign-Ex will attempt to install the required Python modules, but the R libraries must be manually installed. If the required Python modules do not install on their own, they will need to be installed manually as well.

Required Python Modules:
- progress
- networkx; version 1.8.1
- intervaltree

Required R Libraries:
- bestNormalize
- BoutrosLab.plotting.general
- tidyverse

# Set-Up
Quick Start Manual: Follow the instructions to get Benign-Ex up and running and process the provided sample dataset.
User Manual: Follow the instructions to modify Benign-Ex and process your own dataset(s).

# Command to run BENIGNEX:

    python CODE/INTERFACE_BENIGNEX.py PARAMETERS/AUTOMATION

# Citation
If you use the software, please cite: TBD

# Future Improvements
- Support for Python3
