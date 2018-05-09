## genome scan pipeline class file
## generates per SNP Fst input for Yeaman's analyses
## by Christian Sailer
## 6 December 2017
## updated/added step3_nullW on 8 March 2018

import os, sys, statistics, subprocess, argparse
import pandas as pd
import numpy as np
from scipy import stats
from natsort import natsorted

class G_NullW():
	def __init__(self, AFs_table=None, args=None, natsorted=natsorted):
		self.AFs_table = None
		self.args = args
		self.natsored = natsorted

	def step1_nullW(self):
		###### STEP 1 AC TABLE ######
		### Subsets all variant AC sites of the contrast to gene-space
		# set arguments and objects
		args = self.args
		if args.o is 'na':
			outputdir = str(args.i+args.coh1+args.coh2+'/perSNP/')
		else:
			outputdir = str(args.o+args.coh1+args.coh2+'/perSNP/')

		contrast = args.coh1+args.coh2
		self.outputdir = outputdir
		self.contrast = contrast
		logdir = args.o+args.coh1+args.coh2+'/process_log/'
		self.logdir = logdir

		print('\nSTEP 1: Obtain and subset WG AC table\n')

		# Obtain concatenated AC table from outlier scan, a product of G1_outliers.py, variant sites only
		ac = []
		for dirName, subdirList, fileList in os.walk(args.i+contrast):
			for file in fileList:
				if file.startswith(contrast) and file.endswith('WG_ACs'+args.suf+'.txt'):
					ac.append(file)
		ac_sorted = natsorted(ac)

		print('\tFound '+str(len(ac))+' input WG-AC tables.')

		# Convert to bed file
		for table in ac_sorted:
			print('\n\tConvert '+str(table)+' to bed format.')
			basename = table.replace('.txt', '')
			awkcmd = open(outputdir+'subset_AC.unix', 'w')
			awkcmd.write('grep -v CHROM '+str(args.i+contrast+'/'+table)+' | ')
			awkcmd.write("""awk '{OFS="\t"; print $1, $2-1, $2, $3, $4, $5, $6}' """)
			awkcmd.write('> '+outputdir+basename+'.bed')
			awkcmd.close()

			# execute in unix
			cmd = (open(outputdir+'subset_AC.unix', 'r'))
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid, 0)[1]

			# remove up temporary files
			os.remove(outputdir+'subset_AC.unix')
			print('\t'+str(table)+' converted!')

		# subset to genespace (or any other intervall file). the Lyrata gene spce are 34051 genes.
		ac_bed = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.startswith(basename) and file.endswith('.bed'):
					ac_bed.append(file)
		ac_bed_sorted = natsorted(ac_bed)
		print(ac_bed_sorted)

		for table in ac_bed_sorted:
			print('\n\tSubset '+table+' with '+args.genes)
			basename = table.replace('_ACs'+args.suf+'.bed', '')
			bedcmd = open(outputdir+'bed_intersect.unix', 'w')
			bedcmd.write('bedtools intersect -a '+outputdir+table+' -b '+args.genes+' -u ') # -u reports each overlap of a in b only once
			bedcmd.write('> '+outputdir+basename+'_refined_AC'+args.suf+'.table') # no gene name in file 
			bedcmd.close()

			# execute in unix
			cmd = (open(outputdir+'bed_intersect.unix', 'r'))
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid, 0)[1]

			print('\t'+file+' subsetted to '+basename+'_refined_AC'+args.suf+'.table')

			# Obtain gene start, end, length, name from interval bedfile of GS2
			bedgene = open(outputdir+'bed_gene.unix', 'w')
			bedgene.write(""" awk '{FS="\t"; OFS="\t"; print $5, $1, $2, $3, $3-$2}' """)
			bedgene.write(args.genes+ '> '+outputdir+'/temp_'+basename+'_gene_interval.txt')
			bedgene.close()

			# execute in unix
			cmd = open(outputdir+'bed_gene.unix', 'r')
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid,0)[1]

			# remove temorary files
			os.remove(outputdir+'bed_intersect.unix')
			os.remove(outputdir+'bed_gene.unix')

			# add header to bed file, necessary for pasting
			with open(outputdir+basename+'_gene_interval.txt', 'w') as outfile:
				header = 'gene\tscaffold\tstart\tend\tgene_length\n'
				outfile.write(header)
				with open(outputdir+'/temp_'+basename+'_gene_interval.txt') as infile:
					for line in infile:
						outfile.write(line)
			os.remove(outputdir+'/temp_'+basename+'_gene_interval.txt')


	def metric_calculation_gene(self, infile, outfile, winexclcount, file_count, win_count):
		args = self.args
		c1_freq = []
		c2_freq = []
		fst = []
		fst_alt = []

		# Calculate per SNP
		data = infile.readline()
		scaffold, start, position, AC_1, AN_1, AC_2, AN_2 = data.split()
		allele_c1_count = int(AC_1)
		allele_c2_count = int(AC_2)
		allele_c1_freq = allele_c1_count/int(AN_1)
		allele_c2_freq = allele_c2_count/int(AN_2)
		allele_count = allele_c1_count + allele_c2_count
		allele_freq = allele_count/(int(AN_1) + int(AN_2))
		c1_freq.append(allele_c1_freq)
		c2_freq.append(allele_c2_freq)
		c1_minus_c2 = allele_c1_freq - allele_c2_freq
		c2_minus_c1 = allele_c2_freq - allele_c1_freq
		absdiff = abs(allele_c1_freq - allele_c2_freq)
		ht = 2*allele_freq*(1-allele_freq)
		h1 = 2*allele_c1_freq*(1-allele_c1_freq)
		h2 = 2*allele_c2_freq*(1-allele_c2_freq)
		h12 = ((h1*int(AN_1))+(h2*int(AN_2)))/(int(AN_1)+int(AN_2))
		allele_fst = abs(ht-h12)/ht
		# alternate heterozygosity based on autoteraploid frequencies 4p^3q+6p^2q^2+4pq^3
		ht_alt = 4*allele_freq**3*(1-allele_freq) + 6*allele_freq**2*(1-allele_freq)**2 + 4*allele_freq*(1-allele_freq)**3
		h1_alt = 4*allele_c1_freq**3*(1-allele_c1_freq) + 6*allele_c1_freq**2*(1-allele_c1_freq)**2 + 4*allele_c1_freq*(1-allele_c1_freq)**3
		h2_alt = 4*allele_c2_freq**3*(1-allele_c2_freq) + 6*allele_c2_freq**2*(1-allele_c2_freq)**2 + 4*allele_c2_freq*(1-allele_c2_freq)**3
		h12_alt = ((h1_alt*int(AN_1))+(h2_alt*int(AN_2)))/(int(AN_1)+int(AN_2))
		fst_alt = (ht_alt-h12_alt)/ht_alt
		dxy = (allele_c1_freq*(1-allele_c2_freq))+(allele_c2_freq*(1-allele_c1_freq))

		# create outfile
		outfile.write(str(scaffold)+'\t'+
					  str(int(start))+'\t'+
					  str(int(position))+'\t'+
					  str(allele_c1_freq)+'\t'+
					  str(allele_c2_freq)+'\t'+
					  str(absdiff)+'\t'+
					  str(c1_minus_c2)+'\t'+
					  str(c2_minus_c1)+'\t'+
					  str(allele_fst)+'\t'+
					  str(fst_alt)+'\t'+
					  str(dxy)+'\n')
		return winexclcount, file_count, win_count


	def step2_NullW(self):
		###### STEP 2 perSNP ######
		# calculate Fst, allele frequency and allele frequncy difference. Input required as:
		# gene, scaffold, position, cohort1_allele_count, cohort1_all_alleles_count, cohort2_allele_count, cohort2_all_allele_count
		args = self.args
		logdir = self.logdir
		if args.o is 'na':
			outputdir = str(args.i+args.coh1+args.coh2+'/perSNP/')
		else:
			outputdir = str(args.o+args.coh1+args.coh2+'/perSNP/')

		contrast = args.coh1+args.coh2
		self.outputdir = outputdir
		self.contrast = contrast

		print('\nSTEP 2: Calculate Fst per SNP\n')
		bname = contrast+'_WG'
		header = 'scaffold\tstart\tposition\t'+args.coh1+'_freq\t'+args.coh2+'_freq\tabsdiff\t'+args.coh1+'_'+args.coh2+'\t'+args.coh2+'_'+args.coh1+'\tFst_2C\tFst_4C\tDxy\n'
		with open(outputdir+contrast+'_snp_metrics_'+bname+args.suf+'.txt','w') as outfile:
			outfile.write(header)

			ACs = []
			for dirName, subdirList, fileList in os.walk(outputdir):
				for file in fileList:
					if file.endswith('refined_AC'+args.suf+'.table'):
						ACs.append(dirName+file)
			ACs_sort = natsorted(ACs)

			print('\tFound '+str(len(ACs_sort))+' files starting with '+contrast+' and ending with refined_AC'+args.suf+'.table\n')

			# Obtain the list of genes from the gene list text file
			genes = []
			count = 0

			with open(args.gentxt, 'r') as in_gene:
				for gene in in_gene:
					genes.append(gene.replace('\n', ''))
					count +=1
				num_genes = int(count)
				print('\t'+str(count)+' genes in list')
			# Obtain the list of genes from the refined AC table, AC-genes
			for AC in ACs_sort:
				count = 0
				with open(AC, 'r') as infile:
					for line in infile:
						#data = infile.readline()
						scaffold, start, end, AC1, AN1, AC2, AN2 = line.split()
						count +=1
					print('\t'+str(count)+' AC lines read')
							
			file_count = 0
			win_count = 0
			winexclcount = 0

			for AC in ACs_sort:
				# basename = AC.replace('_refined_AC'+args.suf+'.table', '')
				if num_genes > 1:
					file_count += 1
					win_count += num_genes
					noSNP_count = 0

					with open(AC, 'r') as infile:
						for i in range(count):
							winexclcount, file_count, win_count = self.metric_calculation_gene(infile, outfile, winexclcount, file_count, win_count)
				else:
					print('No genes in file.')
		outfile.close()


		print('\n\tAnalyzed '+str(win_count)+' genes')#'+str(file_count)+' files')
		print('\t!! '+str(noSNP_count)+' genes had 0 SNPs!!\n')


	def step3_NullW(self):
		###### STEP 3 SUBSET TO CANDIDATE GENES ######
		args = self.args
		contrast = args.coh1+args.coh2
		bname = contrast+'_WG'
		if args.o is 'na':
			outputdir = str(args.i+args.coh1+args.coh2+'/perSNP/')
		else:
			outputdir = str(args.o+args.coh1+args.coh2+'/perSNP/')

		logdir = outputdir.replace('perSNP', 'process_log/')
		logfile = open(logdir+'PerSNP_'+args.outlier+'.txt', 'w')
		self.logfile = logfile

		print('\nSTEP 3: Restrict SNPs to outlier scan candidate genes\n')

		# create inlist
		# generate bed format file for the next step
		inlist = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.endswith('_snp_metrics_'+bname+args.suf+'.txt'):
					inlist.append(dirName+file)
		in_sort = natsorted(inlist)

		with open(outputdir+contrast+'_snp_metrics_'+bname+args.suf+'.bed','w') as outfile:
			# remove header to turn into bed file
			for file in in_sort:
				count = 0
				with open(file, 'r') as infile:
					infile.readline() # to read and thereby remove the header
					for line in infile:
						outfile.write(line)
						count += 1
					num_snps = count
			logfile.write('Selecting from '+str(num_snps)+' sites\n')
		outfile.close()
		print('\tSelecting from '+str(num_snps)+' sites')
		
		infile = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.endswith('_snp_metrics_'+bname+args.suf+'.bed'):
					infile.append(dirName+file)
		infile = natsorted(infile)

		# subset with refined candidate gene list
		if args.o == 'na':
			outdir = outputdir.replace('perSNP', 'refined')
		else:
			outdir = args.o+contrast+'/refined/'

		for table in infile:
			print('\n\tSubset '+table+' with '+args.outlier+'/'+contrast+'_kept_genes.bed')
			bedcmd = open(outputdir+'bed_intersect.unix', 'w')
			bedcmd.write('bedtools intersect -a '+table+' -b '+outdir+args.outlier+'/'+contrast+'_kept_genes.bed')
			bedcmd.write('> '+outputdir+bname+'_'+args.outlier+'_OS_SNPs_temp'+args.suf+'.txt')  
			bedcmd.close()

			# execute in unix
			cmd = (open(outputdir+'bed_intersect.unix', 'r'))
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid, 0)[1]

		os.remove(outputdir+'bed_intersect.unix')

		with open(outputdir+bname+'_'+args.outlier+'_OS_SNPs'+args.suf+'.bed', 'w') as outfile:
			count_sel = 0
			with open(outputdir+bname+'_'+args.outlier+'_OS_SNPs_temp'+args.suf+'.txt', 'r') as infile:
				for line in infile:
					data = line.split(sep='\t')
					outfile.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\n')
					count_sel += 1
				num_snps2 = count_sel
			logfile.write('Selected '+str(num_snps2)+' outlier SNPs.\n')
		infile.close()
		outfile.close()
		
		# read and write file again to add header
		header = 'scaffold\tstart\tposition\t'+args.coh1+'_freq\t'+args.coh2+'_freq\tabsdiff\t'+args.coh1+'_'+args.coh2+'\t'+args.coh2+'_'+args.coh1+'\tFst_2C\tFst_4C\tDxy\n'
		with open(outputdir+bname+'_'+args.outlier+'_OS_SNPs'+args.suf+'.txt', 'w') as outfile:
			outfile.write(header)
			with open(outputdir+bname+'_'+args.outlier+'_OS_SNPs_temp'+args.suf+'.txt', 'r') as infile:
				for line in infile:
					outfile.write(line)
		infile.close()
		outfile.close()

		# inherit
		self.outputdir = outputdir


	def step4_NullW(self):
		###### STEP 4 MERGE PER SNP METRIC WITH SNPEFF FILE ######
		args = self.args
		logfile = self.logfile
		outputdir = self.outputdir

		print('\nSTEP 4: Add SNPeff annotation to SNPs of outlier genes\n')

		# obtain input files
		vcf_in = []
		vcf_in_path = []
		genes_in = []
		genes_in_path = []

		# get annotated vcf
		for dirName, subdirList, fileList in os.walk('Data/cat_vcfs/'):
			for file in fileList:
				if file.startswith(args.coh1) and file.endswith('_ann'+args.suf+'.txt'):
					vcf_in_path.append(dirName+file)
					vcf_in.append(file)

		#  get SNP intervals from previous step3_NullW()
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.endswith('_OS_SNPs_temp'+args.suf+'.txt'):
					genes_in_path.append(dirName+file)
					genes_in.append(file)

		# turn the annotated text file into bedformat to subset with bedtools before merging
		for vcf in vcf_in_path:
			with open(vcf, 'r') as infile:
				infile.readline() # to remove the header
				with open(outputdir+'temp_ann'+args.suf+'.txt', 'w') as outfile:
					for line in infile:
						ann = line.split(sep="\t")
						outfile.write(ann[0]+'\t'+str(int(ann[1])-1)+'\t'+ann[1]+'\t'+ann[2]+'\t'+ann[3]+'\t'+ann[4]+'\t'+
							ann[5]+'\t'+ann[6]+'\t'+ann[7]+'\t'+ann[8]+'\t'+ann[9]+'\t'+ann[10]+'\t'+ann[11]+'\t'+
							ann[12]+'\t'+ann[13]+'\n')
				outfile.close()
			infile.close()

			for table in genes_in:
				print('\tSubset '+vcf+' with '+table)
				basename = table.replace('_SNPs_temp'+args.suf+'.txt', '')
				bedcmd = open(outputdir+'bed_intersect.unix', 'w')
				bedcmd.write('bedtools intersect -a '+outputdir+table+' -b '+outputdir+'temp_ann'+args.suf+'.txt'+' -wb | ')
				bedcmd.write(""" awk '{FS="\t"; OFS= "\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' """)
				bedcmd.write('> '+outputdir+basename+'_SNPeff_temp'+args.suf+'.txt') # no gene name in file 
				bedcmd.close()

				# execute in unix
				cmd = (open(outputdir+'bed_intersect.unix', 'r'))
				p = subprocess.Popen(cmd, shell=True)
				sts = os.waitpid(p.pid, 0)[1]

				print('\t'+vcf+' subsetted to '+basename+'_SNPeff_temp'+args.suf+'.txt')

		# os.remove(outputdir+'bed_intersect.unix')

		with open(outputdir+basename+'_SNPeff_temp'+args.suf+'.txt', 'r') as infile:
			with open(outputdir+basename+'_SNPeff'+args.suf+'.txt', 'w') as outfile:
				outfile.write('chrom\tstart\tend\t'+args.coh1+'_freq\t'+args.coh2+'_freq\tabsdiff\t'+
				args.coh1+'_'+args.coh2+'\t'+args.coh2+'_'+args.coh1+'\tFst_2C\tFst_4C\tDxy\tRef\tAlt\t'+
				'effect\timpact\tgene\tgene_id\tbiotype\tHGVS_C\tHGVS_D\tcDNA_pos\tCDS_pos\tprotein_pos\n')
				for line in infile:
					outfile.write(line)
			outfile.close()
		infile.close()

		# os.remove(outputdir+'temp_ann'+args.suf+'.txt')

	def step1to4(self):
		self.step1_nullW()
		self.step2_NullW()
		self.step3_NullW()
		self.step4_NullW()

	def step3to4(self):
		self.step3_NullW()
		self.step4_NullW()


if __name__ == '__main__':
	import os, sys