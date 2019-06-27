## genome scan pipeline class file
## has to be loaded within the Gx_x.py genome scan scripts
## by Christian Sailer with help from Katie Barr, 5 December 2016
## updated 24 October 2017, fixed pi and Tajima's D calculation
## add removal of all non-variant sites, 2 March 2018
## Removed DP column from input data and from this script, keep non-variant sites, and corrected pi and Watterson's theta, 6 June 2018
## Added flexible input table format

import os, sys, subprocess, statistics, argparse
from natsort import natsorted
import pandas as pd
import numpy as np
from scipy import stats


## create Class 
class Ginterval():
	def __init__(self, AFs_table=None, three_metrics=None, args=None, natsorted=natsorted, outputdir=None):
		self.AFs_table, self.three_metrics = None, None
		self.args = args
		self.natsorted = natsorted
		self.outputdir = outputdir


	def step1_inputfiles(self):
		###### find input tables and extract population names
		print('\nSTEP 1: Find input files\n')
		print('\tSearching input directory for *_raw.table')

		args = self.args
		tintable = []
		tinfiles = []
		for dirName, subdirList, fileList in os.walk(args.i, topdown = False):
			for table in fileList:
				if table.endswith('_raw.table'):
					tintable.append(dirName+'/'+table) # absolute path to file
					tinfiles.append(table) # file name only
		intable = natsorted(tintable)
		infiles = natsorted(tinfiles)

		print('\tFound '+str(len(intable))+' input tables')

		# inherit for downstream steps
		self.infiles = infiles


	def step1a_subset(self):
		###### Convert input AC tables into bed files and subset to interval
		args = self.args
		infiles = self.infiles
		print('\n\n\t\tINTERVAL BED FILE SPECIFIED – ALL ANALYSIS ON SPECIFIED INTERVALS ONLY\n')
		print('\nSTEP 2: Subset raw input data to defined regionsn\n')
		print('\tObtain and subset AC table(s) to candidate gene intervals')

		# obtain interval naming
		tspace = args.interval.split(sep='/')
		space = tspace[len(tspace)-1]
		# define outputdir
		if args.o is 'na':
			outputdir = str('PerPopMetrics_Result_'+space+'/')
		else:
			outputdir = str(args.o+'PerPopMetrics_Result_'+space+'/')
		# check if output directory exists, create it if necessary
		if os.path.exists(outputdir) == False:
			os.mkdir(outputdir)
			print('\n\tCreated directory '+outputdir)

		# Convert to bed file
		for table in infiles:
			print('\t\tConvert '+str(table)+' to bed format.')
			# extract the population name (start of file name, separated by underlines)
			tin = table.split(sep='/')
			inname = tin[len(tin)-1].split(sep='_')[0]
			
			# read in header of file
			with open(args.i+table, 'r') as infile:
				header = infile.readline()
				grepl = header.split(sep='\t')[0]
			infile.close()

			basename = table.replace('.table', '')
			awkcmd = open(outputdir+'subset_AC.unix', 'w')
			awkcmd.write('grep -v '+grepl+' '+str(args.i+table)+' | ')
			awkcmd.write("""awk '{OFS="\t"; print $1, $2-1, $2, $3, $4}' """)
			awkcmd.write('> '+outputdir+'/'+basename+'.bed')
			awkcmd.close()

			# execute in unix
			cmd = (open(outputdir+'subset_AC.unix', 'r'))
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid, 0)[1]

			# remove up temporary files
			os.remove(outputdir+'subset_AC.unix')
			print('\t\t'+str(table)+' converted!')

			intable = basename+'.bed'

			print('\t\tSubset '+intable+' with '+args.interval)

			bedcmd = open(outputdir+'/bed_intersect.unix', 'w')
			bedcmd.write('bedtools intersect -a '+outputdir+intable+' -b '+args.interval+' -wb | ') # -u reports each overlap of a in b only once
			bedcmd.write(""" awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$4,$5,$9}' """)
			bedcmd.write('> '+outputdir+'/'+basename+'_subsetAC'+args.suf+'.table') # no gene name in file 
			bedcmd.close()

			# execute in unix
			cmd = (open(outputdir+'bed_intersect.unix', 'r'))
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid, 0)[1]

			print('\t\t'+table+' subsetted to '+basename+'_subsetAC'+args.suf+'.table\n')

			# remove temorary files
			os.remove(outputdir+'bed_intersect.unix')
			# os.remove(outputdir+'bed_gene.unix')

			# add header to bed file, necessary for pasting
			with open(outputdir+'/'+basename+'_subsetAC'+args.suf+'.txt', 'w') as outfile:
				header = 'scaffold\tstart\tposition\tAC\tAN\n'
				outfile.write(header)
				with open(outputdir+'/'+basename+'_subsetAC'+args.suf+'.table') as infile:
					for line in infile:
						outfile.write(line)

			# inherit for next step
			self.outputdir = outputdir


	def step2_count(self):
		###### STEP 2 ######
		args = self.args

		# define which outputdir
		if args.interval is 'na':
			if args.o is 'na':
				outputdir = str('PerPopMetrics_Result/')
			else:
				outputdir = str(args.o+'PerPopMetrics_Result/')
					# check if output directory exists, create it if necessary
			if os.path.exists(outputdir) == False:
				os.mkdir(outputdir)
				print('\n\tCreated directory '+outputdir)
		else:
			# obtain interval naming
			tspace = args.interval.split(sep='/')
			space = tspace[len(tspace)-1]
			# define outputdir
			if args.o is 'na':
				outputdir = str('PerPopMetrics_Result_'+space+'/')
			else:
				outputdir = str(args.o+'PerPopMetrics_Result_'+space+'/')

		tintable = []
		tinfiles = []

		if args.interval is 'na':
			print('\nSTEP 2: Count sites that are fixed, not segregating and segragting\n')
			print('\tSearching input directory for *_raw.table')
			for dirName, subdirList, fileList in os.walk(args.i, topdown = False):
				for table in fileList:
					if table.endswith('_raw.table'):
						tintable.append(dirName+'/'+table) # absolute path to file
						tinfiles.append(table) # file name only
			intable = natsorted(tintable)
			infiles = natsorted(tinfiles)

		else:
			print('\nSTEP 3: Count sites that are fixed, not segregating and segragting\n')
			for dirName, subdirList, fileList in os.walk(outputdir, topdown = False):
				for table in fileList:
					if table.endswith('_subsetAC'+args.suf+'.txt'):
						tintable.append(dirName+table) # absolute path to file
						tinfiles.append(table) # file name only
			intable = natsorted(tintable)
			infiles = natsorted(tinfiles)


		print('\tFound '+str(len(intable))+' input tables')

		# Write process log file
		print('\tGathering allele count info:\n')
		if os.path.exists(outputdir+'/process_log/') == False:
			os.mkdir(outputdir+'/process_log/')
		logdir = outputdir+'/process_log/'
		self.logdir = logdir

		# count number of different sites
		for i in range(len(intable)):
			# extract the population name (start of file name, separated by underlines)
			tin = intable[i].split(sep='/')
			inname = tin[len(tin)-1].split(sep='_')[0]

			with open(intable[i], 'r') as infile:
				# print(intable[i])
				infile.readline()
				tANmax_1 = []
				if args.interval is 'na':
					for line in infile:
						scaffold, position, AC_1, AN_1 = line.split()
						tANmax_1.append(int(AN_1))
				else:
					for line in infile:
						scaffold, start, position, AC_1, AN_1, gene = line.split()
						tANmax_1.append(int(AN_1))
				ANmax_1 = max(tANmax_1)
				print('\tFor '+inname+': ')
				print('\tANmax = '+str(ANmax_1)+' in cohort '+inname)

				# count site types
				infile.seek(0) # going back to the top of the input file
				with open(outputdir+'/'+inname+'temp_AC'+args.suf+'.table','w') as outfile:
					header = infile.readline() # to read the first line of the file and thereby remove it from the rest of the loop

					count = 0
					target = 100000
					Zero_count = 0
					All_count = 0
					non_variant_count = 0

					if args.interval is 'na':
						outfile.write('CHROM\tPOS\tAC\tAN\tvariant\n')
					for line in infile:
						if args.interval is 'na':
							scaffold, position, AC, AN = line.split()
							if int(AN_1)>=int(round((1-args.per)*ANmax_1,0)):
								total_allele_count = int(AC)
								total_alleles_sampled = int(AN)
								# Count types of variants and add whether they are fixed (NV), non-segrgating (NV) and segregating (V)
								if total_allele_count == total_alleles_sampled: # fixed site
									All_count +=1
									outfile.write(scaffold+'\t'+position+'\t'+AC+'\t'+AN+'\tNV\n')
								elif total_allele_count == 0: # non segegrating site
									Zero_count +=1
									outfile.write(scaffold+'\t'+position+'\t'+AC+'\t'+AN+'\tNV\n')
								else:
									outfile.write(scaffold+'\t'+position+'\t'+AC+'\t'+AN+'\tV\n')
								count += 1
						
						else:
							scaffold, start, position, AC, AN, gene = line.split()
							if int(AN_1)>=int(round((1-args.per)*ANmax_1,0)):
								total_allele_count = int(AC)
								total_alleles_sampled = int(AN)
								# Count types of variants and add whether they are fixed (NV), non-segrgating (NV) and segregating (V)
								if total_allele_count == total_alleles_sampled: # fixed site
									All_count +=1
									outfile.write(scaffold+'\t'+position+'\t'+AC+'\t'+AN+'\tNV\t'+gene+'\n')
								elif total_allele_count == 0: # non segegrating site
									Zero_count +=1
									outfile.write(scaffold+'\t'+position+'\t'+AC+'\t'+AN+'\tNV\t'+gene+'\n')
								else:
									outfile.write(scaffold+'\t'+position+'\t'+AC+'\t'+AN+'\tV\t'+gene+'\n')
								count += 1

					print('\tTotal sites passing filter in file:\t'+str(count))
					print('\tSites with total allele count=0:\t'+str(Zero_count))
					print('\tVariants fixed in both cohorts:\t\t'+str(All_count))
					print('\tVariant/Segregating sites:\t\t'+str(count-Zero_count-All_count)+'\n')

					# write log file with counts
					with open(logdir+inname+'_inputfilter.csv', 'w') as logfile:
						logfile.write('ANmax in cohort '+inname+','+str(ANmax_1)+'\n')
						logfile.write('Total sites passing in file,'+str(count)+'\n')
						logfile.write('Sites with total allele count = 0,\t'+str(Zero_count)+'\n')
						logfile.write('Fixed sites,'+str(All_count)+'\n')
						logfile.write('Segregating sites,'+str(count-Zero_count-All_count))

				if args.interval is not 'na':
					# sort by gene/interval name in unix
					sortcmd = open(outputdir+'/sortAC.unix', 'w')
					sortcmd.write('sort -k 6 '+outputdir+'/'+inname+'temp_AC'+args.suf+'.table > '+outputdir+'/'+inname+'_AC'+args.suf+'.table')
					sortcmd.close()

					# execute in unix
					cmd = (open(outputdir+'/sortAC.unix', 'r'))
					p = subprocess.Popen(cmd, shell=True)
					sts = os.waitpid(p.pid, 0)[1]

					# remove commandline file and temporary file
					os.remove(outputdir+'/sortAC.unix')
					os.remove(outputdir+'/'+inname+'temp_AC'+args.suf+'.table')
			
		# inherit objects
		self.inname = inname


	def metric_calculation(self, infile, outfile, winexclcount, file_count, win_count, num_sites):
		args = self.args
		scaf = []
		pos = []
		pop_freq = []
		theta_pi = []
		theta_h = []
		temp1 = []
		temp2 = []
		num_snps = 0
		site_count = 0
		can_gene = []
		# calculate for each site/SNP
		for i in range(num_sites):
			data = infile.readline()
			if args.interval is 'na':
				scaffold, position, AC, AN_1, variant = data.split()
			else:
				scaffold, position, AC, AN_1, variant, gene = data.split()
				can_gene.append(gene)
			scaf.append(scaffold)
			pos.append(int(position))
			AN = int(AN_1)
			if int(AC) > 0 and int(AC) != int(AN):
				num_snps += 1
			allele_pop_count = int(AC)
			allele_pop_freq = allele_pop_count/int(AN)
			pop_freq.append(allele_pop_freq)
			# Calcuate theta pi based on heterozygosity
			theta_pi.append((2*allele_pop_count*(int(AN)-allele_pop_count))/(int(AN)*(int(AN)-1)))    
			theta_h.append((2*allele_pop_count*allele_pop_count)/(int(AN)*(int(AN)-1)))
			# harmonic number for Watterson's theta
			for j in range((int(AN)-1)):
				temp1.append(1/(j+1))
			sum_rec_1 = sum(temp1)/(i+1) # correct for the number of positions we calcluated this for, '+1', because the index is 0-based

			# to estimate the variance of pi, watterson's theta and Tajima's D, we have to calculate a2, b1, b2, c1, c2
			# harmonic number a2
			for j in range((int(AN)-1)):
				temp2.append(1/((j+1)*(j+1)))
			sum_rec_2 = sum(temp2)/(i+1)
			# b1 & b2
			temp_b1 = (AN+1)/(3*(AN-1))
			temp_b2 = (2*((AN*AN)+AN+3))/(9*AN*(AN-1))
			# c1 & c2
			temp_c1 = temp_b1-(1/sum_rec_1)
			temp_c2 = temp_b2-(AN+2)/(sum_rec_1*AN)+sum_rec_2/(sum_rec_1*sum_rec_1)
		# calculate the means per window
		wstart = min(pos)
		wend = max(pos)
		wlength = wend-wstart
		wmid = wstart+(wlength/2)
		mean_pop_freq = statistics.mean(pop_freq) # allele frequency of alternate base
		# paramters for variance estimatiion
		a1 = sum_rec_1
		a2 = sum_rec_2
		b1 = temp_b1
		b2 = temp_b2
		c1 = temp_c1
		c2 = temp_c2
		e1 = c1/a1
		e2 = c2/(a1*a1+a2)

		# summary statistics
		sum_theta_pi = sum(theta_pi)
		pi = sum(theta_pi)/len(pos) # nucleotide diversity of the window, probability of different base
		theta_w = num_snps/a1
		# sum_theta_h1 = sum(theta_h1)

		# if there is no segregating site within the specified intervals, set Tajima's D to 0
		try:
			tajD = (sum_theta_pi-theta_w)/((e1*num_snps+e2*num_snps*(num_snps-1))**(1/2))
		except ZeroDivisionError:
			tajD = 0

		# write outfile
		if args.interval is 'na':
			outfile.write(scaf[0]+'\t'+
						str(wstart)+'\t'+
						str(wend)+'\t'+
						str(wmid)+'\t'+
						str(wlength)+'\t'+
						str(len(pos))+'\t'+
						str(num_snps)+'\t'+
						# str(args.snps)+'\t'+
						str(mean_pop_freq)+'\t'+
						str(pi)+'\t'+
						str(sum_theta_pi)+'\t'+
						str(theta_w)+'\t'+
						str(tajD)+'\n')
						# str(tajD_False)+'\n')
						# str(sum_theta_h1)+'\n')
		else:
			outfile.write(scaf[0]+'\t'+
						str(wstart)+'\t'+
						str(wend)+'\t'+
						str(wmid)+'\t'+
						str(wlength)+'\t'+
						str(len(pos))+'\t'+
						str(num_snps)+'\t'+
						# str(args.snps)+'\t'+
						str(mean_pop_freq)+'\t'+
						str(pi)+'\t'+
						str(sum_theta_pi)+'\t'+
						str(theta_w)+'\t'+
						str(tajD)+'\t'+
						str(can_gene[0])+'\n')
		return winexclcount, file_count, win_count, num_sites


	def step3_summaryStats(self):
		###### STEP 3 ######
		# calculate allele frequency, pi, Watterson's theta and Tajima's D. Input required as:
		# scaffold, position, cohort1_allele_count, cohort1_all_alleles_count, variant status
		args = self.args
		logdir = self.logdir
		# define outputdir
		if args.o is 'na':
			outputdir = str('PerPopMetrics_Result/')
		else:
			outputdir = str(args.o+'PerPopMetrics_Result/')
		# check if output directory exists, create it if necessary
		if os.path.exists(outputdir) == False:
			os.mkdir(outputdir)
			print('\n\tCreated directory '+outputdir)

		print('\nSTEP 3: Calculate pi, Wattersons theta, Tajimas D')

		# obtain input files	
		AFs = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.endswith('_AC'+args.suf+'.table') == True:
					AFs.append(dirName+'/'+file)
		AFs_sort = natsorted(AFs)

		print('\n\tFound '+str(len(AFs_sort))+' tables ending with AC'+args.suf+'.table\n')
		
		file_count = 0
		for i in range(len(AFs_sort)):
			# extract the population name (start of file name, separated by underlines)
			tin = AFs_sort[i].split(sep='/')
			inname = tin[len(tin)-1].split(sep='_')[0]
			header = 'scaffold\tstart\tend\tmidpoint\tlength_bp\tbp_sequenced\tS_'+inname+'\t'+inname+'_freq\t'+inname+'_pi\t'+inname+'_thetaPi\t'+inname+'_thetaw\t'+inname+'_TajimasD\n'
			with open(outputdir+'/'+inname+'_metrics_'+str(args.snps)+'SNPs'+args.suf+'.txt','w') as outfile:
				outfile.write(header)

				win_count = 0
				winexclcount = 0

				count = 0 
				var_count = 0
				with open(AFs_sort[i],'r') as infile:
					infile.readline() # to remove header from count
					for line in infile:
						count += 1 # counts the total number sites in file
						# I cannot inherit the number from step 2 when changing args.snps. Redo the variant counting to get args.snps == num_var_sites per window
						scaffold, position, AC_1, AN_1, variant = line.split()
						total_allele_count = int(AC_1)
						total_alleles_sampled = int(AN_1)
						if total_allele_count == total_alleles_sampled:
							pass
						elif total_allele_count == 0:
							pass
						else:
							var_count += 1 # counts only segregating sites
					print('\tFor '+inname+':')
					print('\t\tFound '+str(count)+' sites in file')
					print('\t\tFound '+str(var_count)+' segregating sites (S) in file')
					num_win_var = int(var_count/args.snps)
					print('\t\tProcessing '+str(num_win_var)+' windows with '+str(args.snps)+' variants in '+AFs_sort[i]+'\n')

					# write log file with counts
					with open(logdir+inname+'_processed_sites.txt', 'w') as logfile:
						logfile.write('Found '+str(count)+' sites in file\n')
						logfile.write('Found '+str(var_count)+' segregating sites (S) in file\n')
						logfile.write('Processed '+str(num_win_var)+' windows with '+str(args.snps)+' variants in '+AFs_sort[i]+'\n')

					if num_win_var > 0:
						file_count += 1
						win_count += num_win_var

                    # Create a list with window lengths that each contain args.snps variant sites.
						with(open(AFs_sort[i], 'r')) as infile:
							infile.readline()
							total_count = 0
							site_count = 0
							target = args.snps+1
							win_len = []
							snp_count = 0
							# set condition of how many SNPs to run until change
							while snp_count <= target and total_count < count:
								snp_count = 1
								site_count = 0
								for line in infile:
									site_count += 1
									total_count += 1
									scaffold, position, AC_1, AN_1, variant = line.split()
									if variant == 'V':
										snp_count += 1
									if snp_count%target == 0:
										# print(str(snp_count)+' SNPs '+str(site_count)+' sites')
										win_len.append(site_count)
										break
							infile.seek(0) # go back to the start
							infile.readline() # read header

							for i in range(num_win_var):
								num_sites = int(win_len[i])
								winexclcount, file_count, win_count, num_sites = self.metric_calculation(infile, outfile, winexclcount, file_count, win_count, num_sites)

							infile.close()

		print('\nAnalyzed '+str(win_count)+' windows in '+str(file_count)+' files')
		# print('Excluded '+str(winexclcount)+' windows longer than 26560bp\n')


	def step3a_sumstats_interval(self):
		###### STEP 3 ######
		# calculate allele frequency, pi, Watterson's theta and Tajima's D. Input required as:
		# scaffold, position, cohort1_allele_count, cohort1_all_alleles_count, variant status
		args = self.args

		# obtain interval naming
		tspace = args.interval.split(sep='/')
		space = tspace[len(tspace)-1]
		# define outputdir
		if args.o is 'na':
			outputdir = str('PerPopMetrics_Result_'+space+'/')
		else:
			outputdir = str(args.o+'PerPopMetrics_Result_'+space+'/')

		logdir = outputdir+'/process_log/'
		self.logdir = logdir
		
		print('\nSTEP 4: Calculate pi, Wattersons theta, Tajimas D')

		# obtain input files	
		AFs = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.endswith('_AC'+args.suf+'.table') == True:
					AFs.append(dirName+'/'+file)
		AFs_sort = natsorted(AFs)

		print('\n\tFound '+str(len(AFs_sort))+' tables ending with AC'+args.suf+'.table\n')
		
		file_count = 0
		win_count = 0
		winexclcount = 0

		# Obtain the list of intervalse
		genes = []
		count = 0
		# count how many intervals
		with open(args.interval, 'r') as in_gene:
			for gene in in_gene:
				tgene = gene.split('\t')
				genes.append(tgene[3].replace('\n', ''))
				count +=1
			num_genes = int(count)
			genes = natsorted(genes)
			print('\t'+str(num_genes)+' genes/intervals in list')
			# logfile.write(str(num_genes)+' genes/intervals in list\n')
						

		for i in range(len(AFs_sort)):
			# extract the population name (start of file name, separated by underlines)
			tin = AFs_sort[i].split(sep='/')
			inname = tin[len(tin)-1].split(sep='_')[0]
			print(inname)
			# logfile = open(logdir+inname+'_processed_sites.txt', 'w')

			with open(AFs_sort[i],'r') as infile:
				infile.readline() # to remove header from count
				AC_genes = []
				for line in infile:
					scaffold, position, AC_1, AN_1, variant, gene = line.split()
					AC_genes.append(gene)
					count += 1 # counts the total number sites in file
				# AC_genes = natsorted(AC_genes)
				print('\t'+str(count)+' AC lines read')

				# with open(logdir+inname+'_processed_sites.txt', 'w') as logfile:
				# 	logfile.write(str(count)+' AC lines read\n')

			if args.interval is 'na':
				header = 'scaffold\tstart\tend\tmidpoint\tlength_bp\tbp_sequenced\tS_'+inname+'\t'+inname+'_freq\t'+inname+'_pi\t'+inname+'_thetaPi\t'+inname+'_thetaw\t'+inname+'_TajimasD\n'
			else:
				header = 'scaffold\tstart\tend\tmidpoint\tlength_bp\tbp_sequenced\tS_'+inname+'\t'+inname+'_freq\t'+inname+'_pi\t'+inname+'_thetaPi\t'+inname+'_thetaw\t'+inname+'_TajimasD\tgene\n'				
			with open(outputdir+'/'+inname+'_metrics_'+str(args.snps)+'SNPs'+args.suf+'.txt','w') as outfile:
				outfile.write(header)

				if num_genes > 1:
					file_count += 1
					win_count += num_genes
					noSNP_count = 0

					with open(AFs_sort[i], 'r') as infile:
						# outfile2 = open(outputdir+contrast+'_'+baseparname+'_NO_SNPS.txt', 'w')
						# outfile2.write('gene\n')
						# outkeep = open(outputdir+contrast+'_'+baseparname+'_keep_genes.txt', 'w')
						# outkeep.write('gene\n')
						for j in range(num_genes):
							# count how many lines (==sites) in AC table for each gene
							num_sites = sum(1 for x in AC_genes if x == genes[j])
							# print(str(num_sites)+'\tsites in gene '+genes[i])
							# logsnp.write(str(num_snps)+'\tSNPs in gene '+genes[i]+'\n')
							if num_sites > 0:
								winexclcount, file_count, win_count, num_sites = self.metric_calculation(infile, outfile, winexclcount, file_count, win_count, num_sites)
							else:
								print('No site in candidate gene '+genes[j])
								noSNP_count +=1


				else:
					print('No genes in file.')
			outfile.close()


			print('\n\tAnalyzed '+str(win_count)+' genes')#'+str(file_count)+' files')
			print('\t!! '+str(noSNP_count)+' genes had 0 SNPs!!\n')




		###############################
		# Section to test wrong calculation of tajD
		inlist = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.endswith(str(args.snps)+'SNPs'+args.suf+'.txt'):
					inlist.append(dirName+'/'+file)
		inlist_sort = natsorted(inlist)


		for file in inlist_sort:
			tin = file.split(sep='/')
			inname = tin[len(tin)-1].split(sep='_')[0]
			snp_file = pd.read_table(file, header=0)

		# calculate Tajima's D
			delta1 = snp_file[inname+'_thetaPi']-snp_file[inname+'_thetaw']
			snp_file['FalseTajD_'+inname] = delta1/statistics.stdev(delta1)

			snp_file.to_csv(outputdir+'/'+inname+'_WG_'+str(args.snps)+'SNPs_3metrics'+args.suf+'.txt',sep="\t", index=False)
		# end of section to test wrong TajD
		###############################

    # def step5(self):
    #     ###### STEP 5 ######
    #     # define outliers using top percentile cut-off
    #     print('\n\tSTEP 5: Select outlier windows')
    #     args = self.args
    #     if args.o is 'na':
    #         outputdir = str(args.coh1+args.coh2)
    #     else:
    #         outputdir = str(args.o+args.coh1+args.coh2)
    #     logdir = outputdir+'/process_log/'
    #     cwd = os.getcwd()
    #     if args.o is 'na':
    #         outputdir = str(args.coh1+args.coh2)
    #     else:
    #         outputdir = str(args.o+args.coh1+args.coh2)

    #     contrast = args.coh1 + args.coh2
    #     self.outputdir = outputdir
    #     self.contrast = contrast
    #     self.cwd = cwd

    #     out = []
    #     for dirName, subdirList, fileList in os.walk(outputdir):
    #         for file in fileList:
    #             if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_3metrics'+args.suf+'.txt') == True:
    #                 out.append(dirName+'/'+file)
    #     out_sort = natsorted(out)

    #             outlier_values = open(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_descriptive_stats'+args.suf+'.txt', 'w')
    #             outlier_values.write('stat\tAbsDiff\tAFD_1_2\tDxy\tFst\tFst_4C\tDD\twin_length\tTajD_'+args.coh1+'\tTajD_'+args.coh2+'\tFaiWusH_'+args.coh1+'\tFaiWusH_'+args.coh2+'\n'+
    #                                  'min\t'+str(round(min(absafd),4))+'\t'+str(round(min(afd_c1_c2),4))+'\t'+str(round(min(dxy),4))+'\t'+str(round(min(fst),4))+'\t'+str(round(min(fst_4C),4))+'\t'+str(round(min(dd),4))+'\t'+str(round(min(wlength),0))+'\t'+str(round(min(tajD_c1),4))+'\t'+str(round(min(tajD_c2),4))+'\t'+str(round(min(faiwuH_c1),4))+'\t'+str(round(min(faiwuH_c2),4))+'\n'+
    #                                  'mean\t'+str(round(mean_abs,4))+'\t'+str(round(mean_afd,4))+'\t'+str(round(mean_Dxy,4))+'\t'+str(round(mean_Fst,4))+'\t'+str(round(mean_Fst_4C,4))+'\t'+str(round(mean_DD,4))+'\t'+str(round(mean_win,0))+'\t'+str(round(mean_tajd1,4))+'\t'+str(round(mean_tajd2,0))+'\t'+str(round(mean_faiwu1,4))+'\t'+str(round(mean_faiwu2,0))+'\n'+
    #                                  'median\t'+str(round(median_abs,4))+'\t'+str(round(median_afd,4))+'\t'+str(round(median_Dxy,4))+'\t'+str(round(median_Fst,4))+'\t'+str(round(median_Fst_4C,4))+'\t'+str(round(median_DD,4))+'\t'+str(round(median_win,0))+'\t'+str(round(median_tajd1,4))+'\t'+str(round(median_tajd2,0))+'\t'+str(round(median_faiwu1,4))+'\t'+str(round(median_faiwu2,0))+'\n'+
    #                                  'sd\t'+str(round(sd_abs,4))+'\t'+str(round(sd_afd,4))+'\t'+str(round(sd_Dxy,4))+'\t'+str(round(sd_Fst,4))+'\t'+str(round(sd_Fst_4C,4))+'\t'+str(round(sd_DD,4))+'\t'+str(round(sd_win,0))+'\t'+str(round(sd_tajd1,4))+'\t'+str(round(sd_tajd2,0))+'\t'+str(round(sd_faiwu1,4))+'\t'+str(round(sd_faiwu2,0))+'\n'+
    #                                  'CV\t'+str(round(CV_abs,4))+'\t'+str(round(CV_afd,4))+'\t'+str(round(CV_Dxy,4))+'\t'+str(round(CV_Fst,4))+'\t'+str(round(CV_Fst_4C,4))+'\t'+str(round(CV_DD,4))+'\t'+str(round(CV_win,0))+'\t'+str(round(CV_tajd1,4))+'\t'+str(round(CV_tajd2,0))+'\t'+str(round(CV_faiwu1,4))+'\t'+str(round(CV_faiwu2,0))+'\n'+
    #                                  'max\t'+str(round(max(absafd),4))+'\t'+str(round(max(afd_c1_c2),4))+'\t'+str(round(max(dxy),4))+'\t'+str(round(max(fst),4))+'\t'+str(round(max(fst_4C),4))+'\t'+str(round(max(dd),4))+'\t'+str(round(max(wlength),0))+'\t'+str(round(max(tajD_c1),4))+'\t'+str(round(max(tajD_c2),0))+'\t'+str(round(max(faiwuH_c1),4))+'\t'+str(round(max(faiwuH_c2),0))+'\n'+
    #                                  'Cutoff\t'+str(round(out_abs_up,4))+'\t'+str(round(out_afd_up,4))+'\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_fst_4C,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\t'+str(round(out_tajd1,4))+'\t'+str(round(out_tajd2,0))+'\t'+str(round(out_faiwu1,4))+'\t'+str(round(out_faiwu2,0))+'\n'+
    #                                  'Cutoff2\t'+str(round(out_abs_low,4))+'\t'+str(round(out_afd_low,4))+'\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_fst_4C,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\t'+str(round(out_tajd1,4))+'\t'+str(round(out_tajd2,0))+'\t'+str(round(out_faiwu1,4))+'\t'+str(round(out_faiwu2,0))+'\n')


	def step4_graphs(self):
		###### STEP 4 ######
		## create histograms
		args = self.args

		# define which outputdir
		if args.interval is 'na':
			if args.o is 'na':
				outputdir = str('PerPopMetrics_Result/')
			else:
				outputdir = str(args.o+'PerPopMetrics_Result/')
			# check if output directory exists, create it if necessary
			if os.path.exists(outputdir) == False:
				os.mkdir(outputdir)
				print('\n\tCreated directory '+outputdir)
		else:
			# obtain interval naming
			tspace = args.interval.split(sep='/')
			space = tspace[len(tspace)-1]
			# define outputdir
			if args.o is 'na':
				outputdir = str('PerPopMetrics_Result_'+space+'/')
			else:
				outputdir = str(args.o+'PerPopMetrics_Result_'+space+'/')

		cwd = os.getcwd()
		if args.interval is 'na':
			print('\nSTEP 4: Create population genetics histograms\n')
		else:
			print('\nSTEP 5: Create population genetics histograms\n')
		graphs = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				# if 'metrics' in filename == True:
				if file.endswith(str(args.snps)+'SNPs'+args.suf+'.txt') == True:
					graphs.append(dirName+'/'+file)
		graphs_sorted = natsorted(graphs)

		if os.path.exists(outputdir+'/graphs') == False:
			os.mkdir(outputdir+'/graphs')
		
		print('\tCreating histograms for pi, Tajimas D and window length of populations:')		
		# for i in range(len(graphs_sorted)):
		# 	tin = graphs_sorted[i].split(sep='/')
		# 	inname = tin[len(tin)-1].split(sep='_')[0]
			# print('\t\t'+str(inname))
		for file in graphs_sorted:
			tin = file.split(sep='/')
			inname = tin[len(tin)-1].split(sep='_')[0]
			print(inname)
			count = 0
			rfile = open(outputdir+'/graphs/'+inname+'_graphs.r', 'w')
			rfile.write('# load table\n'+
						'test <- read.delim("'+cwd+'/'+file+'")\n'+
						'#str(test)\n'+
						'#load library\n'+
						'library(ggplot2)\n'+
						'library(methods)\n'+
						'library(grid)\n'+
						'test$pos <- paste(test$scaffold, test$midpoint, sep="-")\n'+
						'p.pi <- qplot('+str(inname)+'_pi, data=test, binwidth=0.001, xlim=c(0,0.5)) + theme_bw()\n'+ #+ geom_vline(xintercept='+str(self.out_afd_up)+', color="grey30", linetype="dashed") + geom_vline(xintercept='+str(self.out_afd_low)+', color="grey30", linetype="dotted") + theme_bw()\n'+
						'p.TajD <- qplot('+str(inname)+'_TajimasD, data=test, binwidth=0.01, xlim=c(-3, +3)) + theme_bw()\n'+ #+ geom_vline(xintercept='+str(self.out_dxy)+', color="grey30", linetype="dashed") + theme_bw()\n'+
						'p.length <- qplot(length_bp, data=test, binwidth=100, xlim=c(0,30000)) + theme_bw()\n\n'+
						'# export as pdf\n'+
						'pdf(file="'+outputdir+'/graphs/'+inname+'_'+str(args.snps)+'SNPs_pi'+args.suf+'.pdf")\n'+
						'p.pi\n'+
						'dev.off()\n'+
						'pdf(file="'+outputdir+'/graphs/'+inname+'_'+str(args.snps)+'SNPs_TajimasD'+args.suf+'.pdf")\n'+
						'p.TajD\n'+
						'dev.off()\n'+
						'pdf(file="'+outputdir+'/graphs/'+inname+'_'+str(args.snps)+'SNPs_length'+args.suf+'.pdf")\n'+
						'p.length\n'+
						'dev.off()\n\n'+
						'# create overview graph\n'+
						'layout <- theme(legend.position="none", axis.text.x=element_blank())\n'+
						'pgenom <- qplot(pos, '+inname+'_pi, data=test, color=scaffold, ylim=c(0,0.5)) + geom_hline(yintercept=median(test$'+inname+'_pi)) + layout\n'+
						'phist <- p.pi + coord_flip() + theme_bw() + geom_vline(xintercept=median(test$'+inname+'_pi))\n'+
						'pdf("'+outputdir+'/graphs/'+inname+'_pi.pdf", width=12, height=4)\n'+
						'grid.newpage()\n'+
						'pushViewport(viewport(layout=grid.layout(1,8)))\n'+
						'vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}\n'+
						'print(pgenom, vp=vplayout(1,1:6))\n'+
						'print(phist, vp=vplayout(1,7:8))\n'+
						'dev.off()\n'+
						'pgenom <- qplot(pos, '+inname+'_TajimasD, data=test, color=scaffold, ylim=c(-3,3)) + geom_hline(yintercept=median(test$'+inname+'_TajimasD)) + layout\n'+
						'phist <- p.TajD + coord_flip() + theme_bw() + geom_vline(xintercept=median(test$'+inname+'_TajimasD))\n'+
						'pdf("'+outputdir+'/graphs/'+inname+'_TajimasD.pdf", width=12, height=4)\n'+
						'grid.newpage()\n'+
						'pushViewport(viewport(layout=grid.layout(1,8)))\n'+
						'vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}\n'+
						'print(pgenom, vp=vplayout(1,1:6))\n'+
						'print(phist, vp=vplayout(1,7:8))\n'+
						'dev.off()\n'+
						'pgenom <- qplot(pos, length_bp, data=test, color=scaffold, ylim=c(0,30000)) + geom_hline(yintercept=median(test$length_bp)) + layout\n'+
						'phist <- p.length + coord_flip() + theme_bw() + geom_vline(xintercept=median(test$length_bp))\n'+
						'pdf("'+outputdir+'/graphs/'+inname+'_length_bp.pdf", width=12, height=4)\n'+
						'grid.newpage()\n'+
						'pushViewport(viewport(layout=grid.layout(1,8)))\n'+
						'vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}\n'+
						'print(pgenom, vp=vplayout(1,1:6))\n'+
						'print(phist, vp=vplayout(1,7:8))\n'+
						'dev.off()\n')
			rfile.close()

			cmd = ('Rscript '+outputdir+'/graphs/'+inname+'_graphs.r')
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid, 0)[1]

			count += 1
			# os.remove(outputdir+'/graphs/'+inname+'_graphs.r')
		print('\n\tDONE\n')    
    








	# define the sequence of steps that should be run
	# def step1to1(self):
	# 	self.step1_inputfiles()
	# def step1to2(self):
	# 	self.step1_inputfiles()
	# 	self.step2_count()
	def step1to4(self):
		self.step1_inputfiles()
		self.step2_count()
		self.step3_summaryStats()
		self.step4_graphs()

	# def step2to4(self):
	# 	self.step2_count()
	# 	self.step3_summaryStats()
	# 	self.step4_graphs()

	def step3to4(self):
		self.step3_summaryStats()
		self.step4_graphs()
	
	def step4to4(self):
		self.step4_graphs()

	def stepsubsetto4(self):
		self.step1_inputfiles()
		self.step1a_subset()
		self.step2_count()
		self.step3a_sumstats_interval()
		self.step4_graphs()

	def stepsubset3to4(self):
		self.step3a_sumstats_interval()
		self.step4_graphs()



	def step2to2(self):
		# self.step1()
		self.step2()

	# def step1to6(self):
	#     self.step1()
	#     self.step2()
	#     self.step3()
	#     self.step4()
	#     self.step5()
	#     self.step6()
	# def step4to6(self):
	#     self.step4()
	#     self.step5()
	#     self.step6()
	# def step5to6(self):
	#     self.step5()
	#     self.step6()

if __name__ == '__main__':
    import os, sys