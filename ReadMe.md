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
scripts are executed when referring to the appropriate script in that folder. See
examples below.
To install python3 packages, on Mac, use homebrew (https://brew.sh) to install python3 via
`brew install python3`. To install any of the modules type e.g. `pip3 install pandas`. If you have several
versions of python you can specify which one you want the module to be installed for,
eg. `pip3.4 install pandas` installs the module 'pandas' for python 3.4.
To install SnpEff go to http://snpeff.sourceforge.net
On Mac you can install it using `brew install snpeff`. However, this might not be the latest version.

## CITATION
If you use these scripts, please cite  
>Preite Veronica, Sailer Christian, Syllwasschy Lara, Bray Sian, Ahmadi Hassan, Krämer Ute and Yant Levi; Convergent evolution in Arabidopsis halleri and Arabidopsis arenosa on calamine metalliferous soils; Phil. Trans. R. Soc. B **374** (1777)

## Preparing input files
After variant calling and filtering using the pipeline of your choice, create per population
tables that contain:  
**CHROM	POS	AC	AN**  
as tab-delimited file, with a filename ending on `_raw.table`. The beginning of the file
has to be your population name, separated from the rest of the file name by `_`.
e.g. `Klet_variants_raw.table`.


## Formulae
Describing populations
We calculate pi based on heterozygosity and calculate Watterson's Theta as described in
'Principles of Population Genetics'. For Tajima's D, we follow Tajima 1989.
```
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

AC...alternate allele count
S....number of segregating sites (called variant, AC)
n....number of sequenced chromosomes (AN)
```
Population differentiation
We calculate Fst based on heterozygosity. For diploid individuals:
```
allele_freq = AC/AN
ht = 2*allele_freq*(1-allele_freq)
h1 = 2*allele_c1_freq*(1-allele_c1_freq)
h2 = 2*allele_c2_freq*(1-allele_c2_freq)
h12 = ((h1*int(AN_1))+(h2*int(AN_2)))/(int(AN_1)+int(AN_2))
allele_fst = abs(ht-h12)/ht
```
Population differentiation
We calculate Fst based on heterozygosity. For tetraploid individuals:
```
ht_alt = 4*allele_freq**3*(1-allele_freq) + 6*allele_freq**2*(1-allele_freq)**2 + 4*allele_freq*(1-allele_freq)**3
h1_alt = 4*allele_c1_freq**3*(1-allele_c1_freq) + 6*allele_c1_freq**2*(1-allele_c1_freq)**2 + 4*allele_c1_freq*(1-allele_c1_freq)**3
h2_alt = 4*allele_c2_freq**3*(1-allele_c2_freq) + 6*allele_c2_freq**2*(1-allele_c2_freq)**2 + 4*allele_c2_freq*(1-allele_c2_freq)**3
h12_alt = ((h1_alt*int(AN_1))+(h2_alt*int(AN_2)))/(int(AN_1)+int(AN_2))
fst_alt = (ht_alt-h12_alt)/ht_alt
```
Dxy and DD:
```
dxy = (allele_c1_freq*(1-allele_c2_freq))+(allele_c2_freq*(1-allele_c1_freq))

absdiff = abs(allele_c1_freq-allele_c2_freq)
prediction = cor_intercept + (cor_slope*absdiff)
DD = pi-prediction
```


# PerPopulationDescription
### SYNTAX
`python3 PerPopMetrics.py -i -snps -per -o -suf`  
If you need help, type  
`python3 PerPopMetrics.py -h`

optional arguments [default value]:  
  `-h`, `--help`           	 	show this help message and exit  
  `-i` inputdir_path     	 	REQUIRED: Full or relative path to input AC table directory  
  `-snps` snps_per_window   Number of variant sites (SNPs) per window [100]  
  `-per` percent_missing_data	Ratio of missing data allowed, eg 0.2 allows 80percent missing data [0.25]  
  `-o` outputdir_path     	Optional output directory relative path, results directory will be a subdir of this one  
  `-suf` suffix_species   	Suffix to append to the file name, e.g. 2-letter species abbreviation, e.g. `_Aa`


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
The DivergenceScan pipeline consists of four parts and takes allele count tables (see *Preparing Input Files*) for two contrasting populations as input. It returns one summary table with Dxy, Fst, DD of all sites and several summary tables that contain only certain outliers and a combination of them (double and triple) only if such were found. Furthermore, it calculates site frequency spectrum metrics (pi, theta). Additionally, it outputs a list with candidate genes, a list with the *A. thaliana* orthologs gene annotation plus the histograms of the popgen metrics (population genetics summary statistics). Moreover, it creates a series of graphs around the candidate loci.  
Works with missing data.  
Most of the steps create logfiles that contain the printed to screen messages as a `.txt`  
file, placed in the folder `process_log`. If run on a SLURM cluster, the printed to screen
messages will be written into the `*.out` files.


### Necessary files
- gene annotation file
- gene ortholog list and gene onthology file


## Description of parts
The four parts are:  
- DS1_outliers.py   Prepares the input AC tables and ends with tables of metrics and outliers, and histograms of metrics. Requires `DS_classes_pi.py` that contains all functions  
- DS2_genes.py 			Finds the genes that overlap with the identified outlier windows  
- DS3_perSNP.py			Produces AFD graphs and popgen metrics graphs for outlier gene lists
- G3_graphs_percentile.py Creates graphical output of results

The three steps have to be run serial as they are dependent on the output of the previous one.

`DS1_outliers.py`, `DS3_perSNP.py` require a class file that has to be in the
same directory as the execution scripts.

### SYNTAX
`python3 DS1_outliers.py -i -coh1 -coh2 -o -snps -cut -per -suf -ploidy`  
`-i`		REQUIRED: Relative or full path to directory that contains the AC tables
		(output of GATKs VariantsToTable) per scaffold or per genome  
`-coh1`	REQUIRED: Name of cohort 1, spelled as in the input tables (literal string, case
		sensitive)  
`-coh2`	REQUIRED: Name of cohort 2, spelled as in the input tables (literal string, case
		sensitive)  
`-o`		Relative path to outputdir, if not specified then the contrast subfolder will be
		created in the current working directory  
`-snps`	Number of SNPs per analysed window, default is 100SNPs  
`-cut`	Top percentile that defines the outlier threshold, default is 0.5  
`-per`	Percentage of allowed missing data specifed as a ratio, eg. for 20% type 0.2,
		default is 0.25  
`-suf`	Suffix to be added at the end, e.g. species shorthand as `_Aa`    
`-ploidy`	For correct heterozygosity estimation of polyploids, this has to be specified. So
		far this has been implemented for autotetraploids. Default is 2.

#### Notes
First, the script tests whether `-o` has been specified and whether the output-subdirectory exists; if not the script creates it. The script runs all steps only the first time for a given setting. Once you change the `-snps` argument, the script does not repeat data preparation and filtering but jumps to step 4. If you change the `-cut` argument after the first run, the script jumps straight to step 5. Histograms are only produced upon changing the `-snps` argument. The outlier combinations are inclusive, that is, e.g. Fst contains all outlier windows that are listed in FstDD.

This first part runs 5 steps within the script:  
step 1: Paste AC tables next to each other to form CHROM1 POS1 AC1 AN1 CHROM2 POS2 AC2 AC2  
step 2: Remove invariant and therefore uninformative sites. The input can contain missing data. Removes sites that have more than the specified percentage of missing data.  
step 3: Compute population genetics summary statistics (popgen metrics) based on non-overlapping SNP-based windows. We do this to avoid inflation of Fst due to a small number of SNPs in a physical window.  
step 4: Determine outliers and select multiple outliers based on cutoff (-cut)  
step 5: Create histograms as pdf of Dxy, Fst, DD, bp length of SNP-based windows

Important output files to look at:  
```
process_log/*
*descriptive_stats.txt
```


`python3 G2_genes.py -i -an -gf -ovlp -suf`  
`-i`		REQUIRED: Relative or full path to outputdir of G1_outliers.py, e.g. `halleri/KromKosi/`  
`-an`		Relative or full path to genome annotation file  
`-gf`		Relative or full path to gene function file of Arabidopsis thaliana orthologs  
`-ovlp`	Percentage of base pairs in the outlier window that should overlap in the genes, default is 0.00001 (equals 1bp overlap)  
`-suf`	Suffix at the end of the filename, has to be identical to the one defined in `DS1_outliers_pi.py`

#### Notes
The script searches the input directory for bedfiles, creates a subdirectory `genes/` and places the gene lists and *Arabidopsis thaliana* gene function lists per outlier combination of popgen metrics there. Additionally it produces one list containing `ALL` in the filename. Here, all candidates genes are listed irrespective of the outlier combination. Use this file to obtain the number of candidate genes.  
Four types of output files are created (example for -ovlp 51):
- `*51ol_genes.gff` 	Subset of the annotation file corresponding to the candidate genes  
- `*51ol_intervals.bed`	Interval list in bed format for candidate genes  
- `*51ol_genes.txt`		Candidate gene list  
- `*51ol_genes_GF.txt`	Gene function of candidate genes of the corresponding *A. thaliana* orthologues


`python3 G3_graphs_percentile.py -coh1 -coh2 -geneor -snps -cut -win -i -ovlp -outlier -suf -ploidy`

`-coh1` REQUIRED: Name of cohort 1, spelled as in the input tables (literal string, case sensitive)  
`-coh2`		REQUIRED: Name of cohort 2, spelled as in the input tables (literal string, case sensitive)  
`-geneor`		Relative of full path to gene orientation file (e.g. Lyrata2genesOriented.out)  
`-snps`		SNP window size in SNPs for which the graphs shall be generated, default is 100  
`-cut`		Top percentile that defines the outlier threshold, default is 0.5  
`-win`		Genome window size of graph in bp, default is 50000  
`-i`			Relative or full path to contrast output directory of part DS1_outliers.py  
`-ovlp`		Percentage of base pairs that should overlap in the search, default is 0.00001, necessary for search pattern  
`-outlier`	Outlier/Outlier combinations for which to produce graphs, possiblities are: Dxy, Fst, DD, DxyFst, DxyDD, FstDD, DxyFstDD  
`-suf`		Suffix to append to file name, eg `_Aa` for *Arabidopsis arenosa*, has to be identical to the one defined in G1_outliers.py, necessary for file search  
`-ploidy`		For correct heterozygosity estimation of polyploids, this has to be specified. So far this has been implemented for autotetraploids. Default is 2.  

CAVEAT: This script runs the main code in R and has to read in all files for every graph generated, which makes it slow.  
The script creates subdirectories within the `genes/` directory according to the outlier combination. For each candidate gene one overview graph with allele frequency difference per SNP, Dxy, Fst, DD and allele frequency difference per SNP window is generated. The gene models are also added as arrows and the candidate gene is highlighted.


# EXAMPLE (incomplete)
Given two cohorts, Klet and Kowa, we have AC tables in the directory `AC_tables`:  
`Klet_raw.table` and `Kowa_raw.table`. Within subdirectories of the `AC_tables` directory:
`AC_tables/scaffold_1/Klet_raw.table`
```
cwd/
	G1_outliers.py
	G2_genes.py
	G3_graphs.py
cwd/AC_tables/
		scaffold_1/
		scaffold_2/
cwd/AC_tables/scaffold_1/
			Klet_raw.table
			Kowa_raw.table
cwd/AC_tables/scaffold_2/
			Klet_raw.table
			Kowa_raw.table
```

In our current working directory (cwd) we have the `AC_tables/` directory and the three python3 scripts. We want to save the DivergenceScan results in the folder `DS_1/`.

`python3 DS1_outliers.py -coh1 Klet -coh2 Kowa -i AC_tables/ -o DS_1/`

This will create the subfolder `DS_1/`, compute the popgen metrics for 100SNP windows and determine the outliers based on the set `-cut` cutoff. Now we have the additional `DS_1/` folder in our cwd. Within it we have the following:
```
cwd/DS_1/
	KletKowa/
cwd/DS_1/KletKowa/
		KletKowa_scaf_1_AFs.table
		KletKowa_scaf_2_AFs.table
		KletKowa_WG_100SNPs_3metrics.txt
		KletKowa_100SNPs_4sd_DxyFstDD_outlieres.csv
		KletKowa_100SNPs_4sd_Fst_outliers.csv
		KletKowa_100SNPs_4sd_DxyFstDD_outliers.bed
		KletKowa_100SNPs_4sd_Fst_outliers.bed
cwd/DS_1/KletKowa/
		graphs/
cwd/DS_1/KletKowa/graphs/
			KletKowa_100SNPs_DD_histogram.pdf
			KletKowa_100SNPs_Dxy_histogram.pdf
			KletKowa_100SNPs_Fst_histogram.pdf
			KletKowa_100SNPs_length_histogram.pdf
```

 You can see that we have `Fst` and `DxyFstDD` outlier lists but no other combinations. The histograms show the distribution and the set cutoff. For the next step, we need the bedfiles and the `annotations/` folder.

```
 cwd/annotations/
 		LyV2.gff
 		LyV2_TAIR10orth_des_20150927.txt
 		Lyrata2GenesOriented.out
```

`python3 DS2_genes.py -i DS_1/KletKowa/ -an annotations/LyV2.gff -gf annotations/LyV2_TAIR10orth_des_20150927.txt`

This creates the subfolder genes within the contrast subfolder:
```
cwd/DS_1/KletKowa/
	graphs/
	lax_genes/
cwd/DS_1/KletKowa/lax_genes/
		KletKowa_100SNPs_4sd_DxyFstDD_0ol_genes.txt
		KletKowa_100SNPs_4sd_DxyFstDD_0ol_genes_GF.txt
		KletKowa_100SNPs_4sd_DxyFstDD_0ol_genes.gff
		KletKowa_100SNPs_4sd_DxyFstDD_0ol_intervals.bed
		KletKowa_100SNPs_4sd_Fst_0ol_genes.txt
		KletKowa_100SNPs_4sd_Fst_0ol_genes_GF.txt
		KletKowa_100SNPs_4sd_Fst_0ol_genes.gff
		KletKowa_100SNPs_4sd_Fst_0ol_intervals.bed
```

Or we request 51% overlap between the SNP windows and the gene annotation:

`python3 G2_genes.py -i GS_1/KletKowa/ -an annotations/LyV2.gff -gf annotations/LyV2_TAIR10orth_des_20150927.txt -ovlp 51`

This creates the subfolder genes within the contrast subfolder:
```
cwd/DS_1/KletKowa/
	graphs/
	genes/
cwd/DS_1/KletKowa/genes/
		KletKowa_100SNPs_4sd_DxyFstDD_51ol_genes.txt
		KletKowa_100SNPs_4sd_DxyFstDD_51ol_genes_GF.txt
		KletKowa_100SNPs_4sd_DxyFstDD_51ol_genes.gff
		KletKowa_100SNPs_4sd_DxyFstDD_51ol_intervals.bed
		KletKowa_100SNPs_4sd_Fst_51ol_genes.txt
		KletKowa_100SNPs_4sd_Fst_51ol_genes_GF.txt
		KletKowa_100SNPs_4sd_Fst_51ol_genes.gff
		KletKowa_100SNPs_4sd_Fst_51ol_intervals.bed
```

The `*_GF.txt` files contain the gene function of the *A. thaliana* ortholog.

As next step, we want to create graphs to annotate the gene list results.

`python3 G3_graphs.py -coh1 Klet -coh2 Kowa -i GS_1/ -geneor annotations/Lyrata2GenesOriented.out -out DxyFstDD -suffix _Aa`

This creates a subfolder `plots_KletKowa_100SNPs_4sd_DxyFstDd_0ol/` within the `genes/` folder. Within this folder, the graphs contain the candidate gene name, the contrast, the SNP window size, the cutoff value, the outlier type and the plotwindowsize, eg:
`AL1G45820_KletKowa_100SNPs_4sd_DxyFstDD_0ol_100kb_Aa`
