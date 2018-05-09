### GENOME SCAN PIPELINE YANT LAB 2016 ###
### individual scripts written by Jeff DaCosta (parts 1-4) streamlined and extended by Christian Sailer 
### 9 January 2017
### updated 7 March 2018

### SHORT DESCRIPTION
This is a command line based tool, developed for use in terminal on Mac.
The genome scan pipeline consists of five parts and takes allele count tables for two 
contrasting populations as input. It returns one summary table with Dxy, Fst, DD of all 
sites and several summary tables that contain only certain outlieres and a combination of 
them (double and triple) only if such were found. Furthermore, it calculates site 
frequency spectrum metrics (pi, theta). Additionally, it outputs a list with candidate 
genes, a list with the A. thaliana orthologs gene annotation plus the histograms of the
popgen metrics (population genetics summary statistics). Moreover, it creates a series of 
graphs around the candidate loci. The last parts are refinements with per gene SNPeff 
annotation and a per SNP metric calculation for selected candidate genes.
Output files are placed into subdirectories named after the contrast.
Works with missing data.
Most of the steps create logfiles that contain the printed to screen messages as a txt 
file, placed in the folder 'process_log'. If run on a SLURM cluster, the printed to screen
messages will be written into the '*.out' files.


### REQUIREMENTS
- gene annotation file
- gene ortholog list and gene onthology file
- python 3
	- numpy
	- pandas
	- argsparse
	- statistics
	- subprocess
	- natsort
	- scipy
	On Mac, use homebrew (https://brew.sh) to install python3 via
	'brew install python3'
	To install any of the modules type eg 'pip3 install pandas'. If you have several 
	versions of python you can specify which one you want the module to be installed for, 
	eg. 'pip3.4 install pandas' installs the module 'pandas' for python 3.4.
- R 3.x
	-ggplot2
- SNPeff
	http://snpeff.sourceforge.net
	On Mac you can install it using
	'brew install snpeff'
	However, this might not be the latest version.
	
Use GATKs VariantsToTable to extract the following fields from cohort vcfs:
CHROM POS AC AN (DP)

The vcfs have to be filtered for:
- Best Practise (BP)
- excess Heterozygosity (HM)
- excess read depth (RDM)

The allele count (AC) table filenames have to contain the cohort name as will be specified
when running the script and end on '_raw.table', e.g. Klet_raw.table for '-coh1 Klet'

To see an example scroll to the end of the document.


### DESCRIPTION OF PARTS
The three parts are:
G1_outliers_pi.py		Prepares the input AC tables and ends with tables of metrics and 
						outliers, and histograms of metrics
						Requires GS_classes_pi.py that contains all functions
G2_genes.py 			Finds the genes that overlap with the identified outlier windows
G3_graphs.py			Produces AFD graphs and popgen metrics graphs for outlier gene 
						lists
G4_refine.py			Calculates metrics per gene and adds the SNPeff annotation for the
						candidate genes
G5_perSNP.py			Calculates metrics per SNP

The five steps have to be run serial as they are dependent on the output of the previous 
one, except G3_graphs.py, which can also be run after G4_refine.py and G5_perSNP.py

G1_outliers.py, G4_refine.py and G5_perSNP.py require a class file that has to be in the
same directory with the scripts.


### G1_outliers_pi.py
#RUN SCRIPT
to see optional arguments and short help message type:
python3 G1_outliers_pi.py -h

syntax:
python3 G1_outliers.py -i -coh1 -coh2 -o -snps -cut -per -suf -ploidy
-i		REQUIRED: Relative or full path to directory that contains the AC tables 
		(output of GATKs VariantsToTable) per scaffold or per genome
-coh1	REQUIRED: Name of cohort 1, spelled as in the input tables (literal string, case 
		sensitive)
-coh2	REQUIRED: Name of cohort 2, spelled as in the input tables (literal string, case 
		sensitive)
-o		Relative path to outputdir, if not specified then the contrast subfolder will be 
		created in the current working directory
-snps	Number of SNPs per analysed window, default is 100SNPs
-cut	Top percentile that defines the outlier threshold, default is 0.5
-per	Percentage of allowed missing data specifed as a ratio, eg. for 20% type 0.2, 
		default is 0.25
-suf	Suffix to be added at the end, e.g. species shorthand as '_Aa'
-ploidy	For correct heterozygosity estimation of polyploids, this has to be specified. So
		far this has been implemented for autotetraploids. Default is 2.

First, the script tests whether -o has been specified and whether the output-subdirectory 
exists; if not the script creates it. The script runs all steps only the first time for a 
given setting. Once you change the -snps argument, the script does not repeat data 
preparation and filtering but jumps to step 4. If you change the -cut argument after the 
first run, the script jumps straight to step 5. Histograms are only produced upon changing
the -snps argument.
The outlier combinations are inclusive, that is, e.g. Fst contains all outlier windows 
that are listed in FstDD.

The first part runs 6 steps within the script:
step 1: Paste AC tables next to each other to form CHROM1 POS1 AC1 AN1 (DP1) CHROM2 POS2 AC2 AC2 (DP2)
step 2: Remove invariant and therefore uninformative sites. The input can contain missing 
		data. Removes sites that have more than the specified percentage of missing data.
(step 3: Remove DP column)
step 4: Compute population genetics summary statistics (popgen metrics) based on 
		non-overlapping SNP-based windows. We do this to avoid inflation of Fst due to a 
		small number of SNPs in a physical window.
step 5: Determine outliers and select multiple outliers based on cutoff (-cut)
step 6: Create histograms as pdf of Dxy, Fst, DD, bp length of SNP-based windows

Importan output files to look at:
process_log/*
*descriptive_stats.txt


### G2_genes.py
#RUN SCRIPT
python3 G2_genes.py -i -an -gf -ovlp -suf

syntax:
-i		REQUIRED: Relative or full path to outputdir of G1_outliers.py, 
		e.g. 'halleri/KromKosi/'
-an		Relative or full path to genome annotation file
-gf		Relative or full path to gene function file of Arabidopsis thaliana orthologs
-ovlp	Percentage of base pairs in the outlier window that should overlap in the genes, 
		default is 0.00001 (equals 1bp overlap)
-suf	Suffix at the end of the filename, has to be identical to the one defined in 
		G1_outliers_pi.py

The script searches the input directory for bedfiles, creates a subdirectory 'genes/' and 
places the gene lists and Arabidopsis thaliana gene function lists per outlier combination
of popgen metrics there. Additionally it produces one list containing 'ALL' in the 
filename. Here, all candidates genes are listed irrespective of the outlier combination. 
Use this file to obtain the number of candidate genes.
Four types of output files are created (example for -ovlp 51):
	*51ol_genes.gff 	Subset of the annotation file corresponding to the candidate genes
	*51ol_intervals.bed	Interval list in bed format for candidate genes
	*51ol_genes.txt		Candidate gene list
	*51ol_genes_GF.txt	Gene function of candidate genes of the corresponding A. thaliana 
						orthologues


### G3_graphs.py
This script can be run after G4_refine.py and G5_perSNP.py
#RUN SCRIPT
python3 G3_graphs.py -coh1 -coh2 -geneor -snps -cut -win -i -ovlp -outlier -suf -ploidy

syntax:
-coh1		REQUIRED: Name of cohort 1, spelled as in the input tables 
			(literal string, case sensitive)
-coh2		REQUIRED: Name of cohort 2, spelled as in the input tables 
			(literal string, case sensitive)
-geneor		Relative of full path to gene orientation file (e.g. Lyrata2genesOriented.out)
-snps		SNP window size in SNPs for which the graphs shall be generated, default is 100
-cut		Top percentile that defines the outlier threshold, default is 0.5
-win		Genome window size of graph in bp, default is 50000
-i			Relative or full path to contrast output directory of part G1_outliers.py
-ovlp		Percentage of base pairs that should overlap in the search, default is 0.00001, 
			necessary for search pattern
-outlier	Outlier/Outlier combinations for which to produce graphs, possiblities are: 
			Dxy, Fst, DD, DxyFst, DxyDD, FstDD, DxyFstDD
-suf		Suffix to append to file name, eg '_Aa' for Arabidopsis arenosa, has to be 
			identical to the one defined in G1_outliers.py, necessary for file search
-ploidy		For correct heterozygosity estimation of polyploids, this has to be specified.
			So far this has been implemented for autotetraploids. Default is 2.

CAVEAT: This script runs the main code in R and has to read in all files for every graph 
generated, which makes it slow.

The script creates subdirectories within the 'genes/' directory according to the outlier 
combination. For each candidate gene one overview graph with allele frequency difference 
per SNP, Dxy, Fst, DD and allele frequency difference per SNP window is generated. The
gene models are also added as arrows and the candidate gene is highlighted.


### G4_refine.py
#RUN SCRIPT
python3 G4_refine.py -i -coh1 -coh2 -snps -o -cut -suf -ovlp -outlier

syntax:
-i			REQUIRED: Full or relative path to directory containing the concatenated AC 
			table (all scaffolds in one table) 
-coh1		REQUIRED: Name of cohort 1, spelled as in the input tables 
			(literal string, case sensitive)
-coh2		REQUIRED: Name of cohort 2, spelled as in the input tables 
			(literal string, case sensitive)
-snps		SNP window size in SNPs for which the graphs shall be generated, 
			default is 100
-cut		Top percentile that defines the outlier threshold, default is 0.5
-o			Relative path to outputdir, if not specified then the subfolder will be 
			created in the current working directory
-ovlp		Percentage of base pairs that should overlap in the search, default is 0.00001, 
			necessary for search pattern
-outlier	Outlier/Outlier combinations for which to generate per SNP metrics, 
			possiblities are: Dxy, Fst, DD, DxyFst, DxyDD, FstDD, DxyFstDD
-suf		Suffix to append to file name, eg '_Aa' for Arabidopsis arenosa, has to be 
			identical to the one defined in G1_outliers.py, necessary for file search

This identifies the gene interval, calculates the popgen metrics array per gene and 
annotates the variant SNPs (filtered in G1_outliers.py) using SNPeff. It then comines the
popgen metrics table with the SNPeff annotation summary per gene.
For SNP-based annotation explore the annotated vcf directly (*ann+suf+.vcf)


### G5_perSNP.py
#RUN SCRIPT
python3 G5_perSNP.py -i -coh1 -coh2 -o -suf -genes -gentxt

syntax:
-i		REQUIRED: Full or relative path to directory containing the concatenated AC 
		table (all scaffolds in one table) 
-coh1	REQUIRED: Name of cohort 1, spelled as in the input tables 
		(literal string, case sensitive)
-coh2	REQUIRED: Name of cohort 2, spelled as in the input tables 
		(literal string, case sensitive)
-o		Relative path to outputdir, if not specified then the subfolder will be 
		created in the current working directory
-suf	Suffix to append to file name, eg '_Aa' for Arabidopsis arenosa, has to be 
		identical to the one defined in G1_outliers.py, necessary for file search
-genes	Full or relative path to genes_only bed file. This bedfile contains the gene 
		intervals.
-gentxt	Full or relative path to all_genes.txt file. This file contains a list of 
		genenames (e.g. AL1G10080), one gene per line.

The output of this step can be used in downstream analyses that require metrics per SNP,
e.g. Sam Yeaman's Null-W.




### EXAMPLE ### (incomplete)
Given two cohorts, Klet and Kowa, we have AC tables in the directory 'AC_tables' and 
Klet_raw.table 
and Kowa_raw.table within subdirectories of the AC_tables directory:
AC_tables/scaffold_1/Klet_raw.table

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

In our current working directory (cwd) we have the 'AC_tables/' directory and the three python3 scripts.
We want to save the genome scan results in the folder 'GS_1/'.

python3 G1_outliers.py -coh1 Klet -coh2 Kowa -i AC_tables/ -o GS_1/

This will create the subfolder 'GS_1/', compute the popgen metrics for 100SNP windows and determine 
the outliers based on mean +/- 4*sd (99.9% CI) cutoff.
Now we have the additional 'GS_1/' folder in our cwd.
Within it we have the following:

cwd/GS_1/
	KletKowa/
cwd/GS_1/KletKowa/
		KletKowa_scaf_1_AFs.table
		KletKowa_scaf_2_AFs.table
		KletKowa_WG_100SNPs_3metrics.txt
		KletKowa_100SNPs_4sd_DxyFstDD_outlieres.csv
		KletKowa_100SNPs_4sd_Fst_outliers.csv
		KletKowa_100SNPs_4sd_DxyFstDD_outliers.bed
		KletKowa_100SNPs_4sd_Fst_outliers.bed
cwd/GS_1/KletKowa/
		graphs/
cwd/GS_1/KletKowa/graphs/
			KletKowa_100SNPs_DD_histogram.pdf
			KletKowa_100SNPs_Dxy_histogram.pdf
			KletKowa_100SNPs_Fst_histogram.pdf
			KletKowa_100SNPs_length_histogram.pdf
 
 You can see that we have 'Fst' and 'DxyFstDD' outlier lists but no other combinations.
 The histograms show the distribution and the 99% and 99.9% CI cutoff.
 For the next step, we need the bedfiles and the 'annotations/' folder.
 
 cwd/annotations/
 		LyV2.gff
 		LyV2_TAIR10orth_des_20150927.txt
 		Lyrata2GenesOriented.out
 		
python3 G2_genes.py -i GS_1/KletKowa/ -an annotations/LyV2.gff -gf annotations/LyV2_TAIR10orth_des_20150927.txt

This creates the subfolder genes within the contrast subfolder:
cwd/GS_1/KletKowa/
	graphs/
	lax_genes/
cwd/GS_1/KletKowa/lax_genes/
		KletKowa_100SNPs_4sd_DxyFstDD_0ol_genes.txt
		KletKowa_100SNPs_4sd_DxyFstDD_0ol_genes_GF.txt
		KletKowa_100SNPs_4sd_DxyFstDD_0ol_genes.gff
		KletKowa_100SNPs_4sd_DxyFstDD_0ol_intervals.bed
		KletKowa_100SNPs_4sd_Fst_0ol_genes.txt
		KletKowa_100SNPs_4sd_Fst_0ol_genes_GF.txt
		KletKowa_100SNPs_4sd_Fst_0ol_genes.gff
		KletKowa_100SNPs_4sd_Fst_0ol_intervals.bed
		
Or we request 51% overlap between the SNP windows and the gene annotation:

python3 G2_genes.py -i GS_1/KletKowa/ -an annotations/LyV2.gff -gf annotations/LyV2_TAIR10orth_des_20150927.txt -ovlp 51

This creates the subfolder genes within the contrast subfolder:
cwd/GS_1/KletKowa/
	graphs/
	genes/
cwd/GS_1/KletKowa/genes/
		KletKowa_100SNPs_4sd_DxyFstDD_51ol_genes.txt
		KletKowa_100SNPs_4sd_DxyFstDD_51ol_genes_GF.txt
		KletKowa_100SNPs_4sd_DxyFstDD_51ol_genes.gff
		KletKowa_100SNPs_4sd_DxyFstDD_51ol_intervals.bed
		KletKowa_100SNPs_4sd_Fst_51ol_genes.txt
		KletKowa_100SNPs_4sd_Fst_51ol_genes_GF.txt
		KletKowa_100SNPs_4sd_Fst_51ol_genes.gff
		KletKowa_100SNPs_4sd_Fst_51ol_intervals.bed


The *_GF.txt files contain the gene function of the A. thaliana ortholog.


As next step, we want to create graphs to annotate the gene list results.

python3 G3_graphs.py -coh1 Klet -coh2 Kowa -i GS_1/ -geneor annotations/Lyrata2GenesOriented.out -out DxyFstDD -suffix _Aa

This creates a subfolder 'plots_KletKowa_100SNPs_4sd_DxyFstDd_0ol/ within the 'genes/' folder
Within this folder, the graphs contain the candidate gene name, the contrast, the SNP window size, the cutoff value, the outlier type and the plotwindowsize, eg:
AL1G45820_KletKowa_100SNPs_4sd_DxyFstDD_0ol_100kb_Aa