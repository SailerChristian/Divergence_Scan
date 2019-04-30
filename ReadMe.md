# DivergenceScan pipeline Yant lab 2018
# Describing populations & DivergenceScan using population genetics summary statistics

## Description
This is a command line based tool, developed for use in terminal on Mac.
The DivergenceScan pipeline consists of two independent parts: 
- Per population description
- DivergenceScan

The PerPopulationDescription is for descriptive purposes only. The scripts divide the provided
input files into windows of a set number of variants and calculates pi and Tajima's D for each
window.
The DivergenceScan script takes the same input files and calculates divergence and diversity
metrics between two cohorts (e.g. two populations). The script divides the data into windows
of a fixed number of SNPs and calculates Fst, dxy and DD and finds the windows that are the 
most diverged by setting an upper percentile cutoff. The script allows for missing data.


## SOFTWARE REQUIREMENTS
- python 3.x
 - package 'argparse'
 - package 'natsort'
 - package 'pandas
 - package 'numpy'
 - package 'scipy'
 - package 'subprocess'

- R 3.x
 - package 'ggplot2'
	
- SNPeff
	
## INSTALLATION
Pull down the code from the repository and save the scripts in a separate folder. The python3
scripts are executed when calling referring to the appropriate script in that folder. See
examples below.
To install python3 packages, on Mac, use homebrew (https://brew.sh) to install python3 via
`brew install python3`. To install any of the modules type eg `pip3 install pandas`. If you have several 
versions of python you can specify which one you want the module to be installed for, 
eg. `pip3.4 install pandas` installs the module 'pandas' for python 3.4.
To install SnpEff go to http://snpeff.sourceforge.net
On Mac you can install it using `brew install snpeff`. However, this might not be the latest version.


## Preparing input files
After variant calling and filtering using the pipeline of your choice, create per population
tables that contain:
**CHROM	POS	AC	AN**
as tab-delimited file, with a filename ending on `_raw.table`. The beginning of the file
has to be your population name, separated from the rest of the file name by `_`.
e.g. `Klet_variants_raw.table`.


## Formulae
We calculate pi based on heterozygosity and calculate Watterson's Theta as described in
'Principles of Population Genetics'. For Tajima's D, we follow Tajima 1989.
```
AC = alternate allele count
Theta(pi) = (2*AC*(n-1))/(n*n-1)
pi = Theta(pi)/len(sequence)
len(sequence) = number of bases passing filters

Theta(w) = S/a1

a1 = sum(1/i) with i=1..n-1
a2 = sum(1/i^2) with i=1..n-1
b1 = (n+1)/(3*(n-1))
b2 = (2*(n^2+n+3))/(9n*(n-1))
c1 = b1-1/a1
c2 = b2-(n+2)/(a1*n)+a2/a1^2
e1 = c1/a1
e2 = c2/(a1^2+a2)

Tajima's D = (Theta(pi)-Theta(w)/sqrt(e1*S+e2*S*(S-1))

S...number of segregating sites (called variant, AC)
n...number of sequenced chromosomes (AN)
```


# PerPopulationDescription
### SYNTAX
`python3 PerPopMetrics.py -i -snps -per -o -suf`
If you need help, type
`python3 PerPopMetrics.py -h`

optional arguments [default value]:
  -h, --help           	 	show this help message and exit
  -i inputdir_path     	 	REQUIRED: Full or relative path to input AC table directory
  -snps snps_per_window 	Number of variant sites (SNPs) per window [100]
  -per percent_missing_data	Ratio of missing data allowed, eg 0.2 allows 80percent missing data [0.25]
  -o outputdir_path     	Optional output directory relative path, results directory will be a subdir of this one
  -suf suffix_species   	Suffix to append to the file name, e.g. 2-letter species abbreviation, e.g. `_Aa`


The output file is a tab-delimited text file that can be used in R, Excel or any other 
software that understands tab-delimited files.
A log directory is created that saves the printed to screen messages.
Additionally, the script produces histograms of pi, Tajima's D and the length of the 
windows, plus an overview graph. The line through the grap is the genome-wide median.


### Script files
`ReadMe_perPopMetric.md`
`PerPopMetrics.py`				Script to execute
`PerPopMetrics_classes_1.py`		This files contains the code for the calculations. Needs to be in the same folder as the executing script


# DivergenceScan