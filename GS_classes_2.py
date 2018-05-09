## genome scan pipeline class file
## for second, hierachial round of genome scanning
## has to loaded with the Gx_x.py genome scan scripts
## by Christian Sailer
## updated 24 October 2017, corrected pi and Tajima's D calculation
## updated 8 March 2018, added logfiles

import os, sys, statistics, subprocess, argparse
import pandas as pd
import numpy as np
from scipy import stats
from natsort import natsorted

class G2():
	def __init__(self, AFs_table=None, args=None, natsorted=natsorted):
		self.AFs_table = None
		self.args = args
		self.natsorted = natsorted
		
	def step1_gene(self):
		###### STEP 1 GENE ######
		# Subset the concatenated AC (allele count) table using the gene list bed files
		args = self.args
		if args.o is 'na':
			outputdir = str(args.coh1+args.coh2+'/refined/'+args.outlier+'/')
		else:
			outputdir = str(args.o+args.coh1+args.coh2+'/refined/'+args.outlier+'/')
		
		contrast = args.coh1 + args.coh2
		logdir = args.o+args.coh1+args.coh2+'/process_log/'
		self.outputdir = outputdir
		self.contrast = contrast
		self.logdir = logdir

		# open logfile
		logfile = open(logdir+'Refine_summary_'+args.outlier+'.txt', 'w')
		print('\n\tSTEP 1 GENE: Convert AC-table to bed file\n')
		## Step 1a: Convert AC table to bed format
		logfile.write('\tSTEP 1 GENE: Convert AC-table to bed file\n')
		# Obtain concatenated filtered AC table from contrast (step G1)
		af = []
		for dirName, subdirList, fileList in os.walk(args.i):
			for file in fileList:
				if file.startswith(contrast) and file.endswith('WG_ACs'+args.suf+'.txt'):
					af.append(file)
		af_sorted = natsorted(af)
		print('Found '+str(len(af))+' input WG-AC table(s).')
		logfile.write('Found '+str(len(af))+' input WG-AC table(s).\n')

		# Obtain gene list bedfile to use for subsetting
		inall = []
		search1 = str('_'+str(args.snps)+'SNPs')
		search2 = str('_'+str(int(10000*args.cut))+'ppm')
		for dirName, subdirList, fileList in os.walk(args.i+'genes/'):
			for file in fileList:
				if search1 in file and search2 in file:
					inall.append(file)
		# restrict to specific outliers
		genes_bed = []
		for file in inall:
			if str('_'+args.outlier+'_'+str(int(args.ovlp))) in file and file.endswith('.bed'):
				genes_bed.append(file)
		genes_bed_sorted = natsorted(genes_bed) # selection intervals
		print('Found '+str(len(genes_bed_sorted))+' outlier interval bed file(s).')
		logfile.write('Found '+str(len(genes_bed_sorted))+' outlier interval bed file(s).\n')
 
		# Conversion
		for table in af_sorted:
			print('\nConvert '+str(table)+' to bed format.')
			basename = table.replace('.txt', '')
			awkcmd = open(outputdir+'subset_AC.unix', 'w')
			awkcmd.write('grep -v CHROM '+str(args.i+'/'+table)+' | ')
			awkcmd.write("""awk '{FS=" "; OFS="\t"; print $1, $2-1, $2, $3, $4, $5, $6}' """)
			awkcmd.write('> '+outputdir+basename+'.bed')
			awkcmd.close()

			# execute in unix
			cmd = (open(outputdir+'subset_AC.unix', 'r'))
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid, 0)[1]

			os.remove(outputdir+'subset_AC.unix')
			print(str(table)+' converted!!')
			logfile.write(str(table)+' converted to bed format.\n')
		# logfile.close()

		self.genes_bed_sorted = genes_bed_sorted
		self.logfile = logfile

	def step2_gene(self):
		###### STEP 2 GENE ######
		# call inherited objects
		args = self.args
		outputdir = self.outputdir
		contrast = self.contrast
		logdir = self.logdir
		logfile = self.logfile

		# Subset the converted AC table and gene interval file
		genes_bed_sorted = self.genes_bed_sorted
		ac_bed = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.startswith(contrast) and file.endswith(args.suf+'.bed'):
					ac_bed.append(file)
		ac_bed_sorted = natsorted(ac_bed)

		# # reopen logfile
		# logfile = open(logdir+'Refine_summary_'+args.outlier+'.txt', 'w')
		print('\n\tSTEP 2 GENE: Subset files\n')
		logfile.write('\n\tSTEP 2 GENE: Subset files\n')
		for table in ac_bed_sorted:
			for file in genes_bed_sorted:
				print('Subset '+str(table)+' with list '+str(file))
				basename = file.replace('_intervals'+args.suf+'.bed', '')
				bedcmd = open(outputdir+'bed_intersect.unix', 'w')
				bedcmd.write('bedtools intersect -a '+args.i+'genes/'+file+' -b '+outputdir+table+' -wb | ')
				bedcmd.write(""" awk '{$1=$2=$3""; print $4,$5,$6,$7,$8,$9,$10,$11}' """)
				bedcmd.write('| ')
				bedcmd.write(""" sed -n -e 's/^.*;Parent=//p' """)
				bedcmd.write('| sort > '+outputdir+basename+'_refined_AC'+args.suf+'.table') # sort list 
				bedcmd.close()

				# execute in unix
				cmd = (open(outputdir+'bed_intersect.unix', 'r'))
				p = subprocess.Popen(cmd, shell=True)
				sts = os.waitpid(p.pid, 0)[1]

				print(file+' subsetted to '+basename+'_refined_AC'+args.suf+'.table')
				logfile.write(file+' subsetted to '+basename+'_refined_AC'+args.suf+'.table\n')

				# Obtain gene start, end, length, name from interval bedfile of GS2
				bedgene = open(outputdir+'bed_gene.unix', 'w')
				bedgene.write(""" sed -n -e 's/.t.*;Parent=.*//p' """)
				bedgene.write(args.i+'genes/'+file+' | sort -u -k 4 | ') # -u sorts uniqu per row, -k4 uses the 4th column to sort
				bedgene.write(""" sed 's/ID=//g' """ '| ')
				bedgene.write(""" awk '{FS="\t"; OFS="\t"; print $4, $1, $2, $3, $3-$2}' """)
				bedgene.write('> '+outputdir+'/temp_'+basename+'_gene_interval.txt')
				bedgene.close()

				# execute in unix
				cmd = open(outputdir+'bed_gene.unix', 'r')
				p = subprocess.Popen(cmd, shell=True)
				sts = os.waitpid(p.pid,0)[1]

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

		# logfile.close()
		# inherit objects
		self.logdir = logdir
		self.logfile = logfile

	
	def metric_calculation_gene(self, infile, outfile, winexclcount, file_count, win_count, num_snps):
		args = self.args
		gene = []
		#scaf = []
		#pos = []
		c1_freq = []
		c2_freq = []
		absdiff = []
		c1_minus_c2 = []
		c2_minus_c1 = []
		theta_pi1 = []
		theta_pi2 = []
		theta_h1 = []
		theta_h2 = []
		temp1 = []
		temp2 = []
		num_snps_1 = 0
		num_snps_2 = 0
		fst = []
		dxy = []
		varu = []
		#print(str(num_spns))
		for i in range(num_snps):
			data = infile.readline()
			can_gene, scaffold, start, position, AC_1, AN_1, AC_2, AN_2 = data.split()
			gene.append(can_gene)
			#scaf.append(scaffold)
			#pos.append(int(position))
			if int(AC_1) > 0:
				num_snps_1 += 1
			if int(AC_2) > 0:
				num_snps_2 += 1
			allele_c1_count = int(AC_1)
			allele_c2_count = int(AC_2)
			allele_c1_freq = allele_c1_count/int(AN_1)
			allele_c2_freq = allele_c2_count/int(AN_2)
			allele_count = allele_c1_count + allele_c2_count
			allele_freq = allele_count/(int(AN_1) + int(AN_2))
			c1_freq.append(allele_c1_freq)
			c2_freq.append(allele_c2_freq)
			c1_minus_c2.append(allele_c1_freq - allele_c2_freq)
			c2_minus_c1.append(allele_c2_freq - allele_c1_freq)
			absdiff.append(abs(allele_c1_freq - allele_c2_freq))
			theta_pi1.append((2*allele_c1_count*(int(AN_1)-allele_c1_count))/(int(AN_1)*(int(AN_1)-1)))    
			theta_pi2.append((2*allele_c2_count*(int(AN_2)-allele_c2_count))/(int(AN_2)*(int(AN_2)-1)))
			theta_h1.append((2*allele_c1_count*allele_c1_count)/(int(AN_1)*(int(AN_1)-1)))
			theta_h2.append((2*allele_c2_count*allele_c2_count)/(int(AN_2)*(int(AN_2)-1)))            
			for j in range((int(AN_1)-1)):
				temp1.append(1/(j+1))
			sum_rec_1 = sum(temp1)
			for k in range((int(AN_2)-1)):
				temp2.append(1/(k+1))
			sum_rec_2 = sum(temp2)
			ht = 2*allele_freq*(1-allele_freq)
			h1 = 2*allele_c1_freq*(1-allele_c1_freq)
			h2 = 2*allele_c2_freq*(1-allele_c2_freq)
			h12 = ((h1*int(AN_1))+(h2*int(AN_2)))/(int(AN_1)+int(AN_2))
			allele_fst = abs(ht-h12)/ht
			fst.append(allele_fst)
			dxy.append((allele_c1_freq*(1-allele_c2_freq))+(allele_c2_freq*(1-allele_c1_freq)))
			varu1 = (1/(int(AN_1)*int(AN_1)))*allele_c1_freq*(1-allele_c1_freq)
			varu2 = (1/(int(AN_2)*int(AN_2)))*allele_c2_freq*(1-allele_c2_freq)
			varu.append(varu1 + varu2)

		#wstart = min(pos)
		#wend = max(pos)
		#wlength = wend-wstart
		#wmid = wstart+(wlength/2)
		mean_c1_freq = statistics.mean(c1_freq)
		mean_c2_freq = statistics.mean(c2_freq)
		mean_absdiff = statistics.mean(absdiff)
		mean_c1_minus_c2 = statistics.mean(c1_minus_c2)
		mean_c2_minus_c1 = statistics.mean(c2_minus_c1)
		# some populations have no snps in a gene, set those to pi=0
		if num_snps_1 > 0:
			sum_theta_pi1 = sum(theta_pi1)/num_snps_1
		else:
			sum_theta_pi1 = sum(theta_pi1)
		if num_snps_2 > 0:
			sum_theta_pi2 = sum(theta_pi2)/num_snps_2
		else:
			sum_theta_pi2 = sum(theta_pi2)
		theta_w1 = num_snps_1/sum_rec_1
		theta_w2 = num_snps_2/sum_rec_2
		sum_theta_h1 = sum(theta_h1)
		sum_theta_h2 = sum(theta_h2)
		mean_fst = statistics.mean(fst)
		sum_dxy = sum(dxy)
		win_dxy = (1/num_snps)*sum_dxy
		mean_varu = statistics.mean(varu)

#             string = "scaf[0]\t %d \t %d \t %d" % (wstart, wend, wmid)
#             outfile.write(string)
		outfile.write(gene[0]+'\t'+
					  #scaf[0]+'\t'+
					  #str(wstart)+'\t'+
					  #str(wend)+'\t'+
					  #str(wlength)+'\t'+
					  str(num_snps)+'\t'+
					  str(mean_c1_freq)+'\t'+
					  str(mean_c2_freq)+'\t'+
					  str(mean_absdiff)+'\t'+
					  str(mean_c1_minus_c2)+'\t'+
					  str(mean_c2_minus_c1)+'\t'+
					  str(mean_varu)+'\t'+
					  str(sum_theta_pi1)+'\t'+
					  str(sum_theta_pi2)+'\t'+
					  str(theta_w1)+'\t'+
					  str(theta_w2)+'\t'+
					  str(sum_theta_h1)+'\t'+
					  str(sum_theta_h2)+'\t'+
					  str(mean_fst)+'\t'+
					  str(win_dxy)+'\n')
		return winexclcount, file_count, win_count, num_snps
			

	def step3_gene(self):
		###### STEP 3 GENE ######
		# calculate Dxy, Fst, pi, allele frequency and allele frequncy difference. Input required as:
		# gene, scaffold, position, cohort1_allele_count, cohort1_all_alleles_count, cohort2_allele_count, cohort2_all_alleles_count
		args = self.args
		logdir = self.logdir
		logfile = self.logfile
		if args.o is 'na':
			outputdir = str(args.coh1+args.coh2+'/refined/'+args.outlier+'/')
		else:
			outputdir = str(args.o+args.coh1+args.coh2+'/refined/'+args.outlier+'/')

		contrast = args.coh1 + args.coh2
		self.outputdir = outputdir
		self.contrast = contrast

		# # reopen logfile
		# logfile = open(logdir+'Refine_summary_'+args.outlier+'.txt', 'w')

		# Obtain gene list text file, only gene names, one per line
		in_gene_all = []
		search1 = str('_'+str(args.snps)+'SNPs')
		search2 = str('_'+str(int(10000*args.cut))+'ppm')
		for dirName, subdirList, fileList in os.walk(args.i+'genes/'):
			for file in fileList:
				if search1 in file and search2 in file:
					in_gene_all.append(file)                    
		# restrict to specific outliers and text file
		gene_list = []
		for file in in_gene_all:
			if str('_'+args.outlier+'_'+str(int(args.ovlp))) in file and file.endswith(args.suf+'.txt'):
				gene_list.append(args.i+'genes/'+file)
		gene_list_sorted = natsorted(gene_list) # to restrict windows

		baseparname = (str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_'+args.outlier+'_'+str(int(args.ovlp))+'ol')
		print('\n\tSTEP 3 GENE: Calculate Dxy, Fst, DD, Pi per population, Tajimas D per pop, thetas')
		logfile.write('\n\tSTEP 3 GENE: Calculate Dxy, Fst, DD, Pi per population, Tajimas D per pop, thetas.\n')
		header = 'can_gene\tnum_snps\t'+args.coh1+'_freq\t'+args.coh2+'_freq\tabsdiff\t'+args.coh1+'_'+args.coh2+'\t'+args.coh2+'_'+args.coh1+'\tvaru\t'+args.coh1+'_pi\t'+args.coh2+'_pi\t'+args.coh1+'_thetaw\t'+args.coh2+'_thetaw\t'+args.coh1+'_thetah\t'+args.coh2+'_thetah\tFst\tDxy\n'
		with open(outputdir+contrast+'_temp_gene_metrics_'+baseparname+args.suf+'.txt','w') as outfile:
			outfile.write(header)

			tACs = []
			ACs = []
			for dirName, subdirList, fileList in os.walk(outputdir):
				for file in fileList:
					if search1 in file and search2 in file:
						tACs.append(file)
			for file in tACs:
				if str('_'+args.outlier+'_'+str(int(args.ovlp))) in file and file.endswith('refined_AC'+args.suf+'.table'):
					ACs.append(dirName+file)
			ACs_sort = natsorted(ACs)

			print('\nFound '+str(len(ACs_sort))+' files starting with '+contrast+' and ending with refined_AC'+args.suf+'.table\n')
			logfile.write('Found '+str(len(ACs_sort))+' files starting with '+contrast+' and ending with refined_AC'+args.suf+'.table\n')

			# Obtain the list of genes from the gene list text file
			for ingene in gene_list_sorted:
				genes = []
				count = 0

				with open(ingene, 'r') as in_gene:
					for gene in in_gene:
						genes.append(gene.replace('\n', ''))
						count +=1
					num_genes = int(count)
					print(str(count)+' genes in list')
					logfile.write(str(count)+' genes in list\n')
			# Obtain the list of genes from the refined AC table, AC-genes
			for AC in ACs_sort:
				AC_genes = []
				count = 0
				with open(AC, 'r') as infile:
					for line in infile:
						#data = infile.readline()
						can_gene, scaffold, start, end, AC1, AN1, AC2, AN2 = line.split()
						AC_genes.append(can_gene)
						count +=1
					print(str(count)+' AC lines read')
					logfile.write(str(count)+' AC lines read\n')
				
			
			file_count = 0
			win_count = 0
			winexclcount = 0

			# create new logfile for SNPs in genes
			logsnp = open(logdir+'SNPs_in_'+args.outlier+'_genes.txt', 'w')
			for AC in ACs_sort:
				basename = AC.replace('_refined_AC'+args.suf+'.table', '')
				if num_genes > 1:
					file_count += 1
					win_count += num_genes
					noSNP_count = 0

					with open(AC, 'r') as infile:
						outfile2 = open(outputdir+contrast+'_'+baseparname+'_NO_SNPS.txt', 'w')
						outfile2.write('gene\n')
						outkeep = open(outputdir+contrast+'_'+baseparname+'_keep_genes.txt', 'w')
						outkeep.write('gene\n')
						for i in range(num_genes):
							# count how many lines (==SNPs) in AC table for each gene
							num_snps = sum(1 for x in AC_genes if x == genes[i])
							print(str(num_snps)+'\tSNPs in gene '+genes[i])
							logsnp.write(str(num_snps)+'\tSNPs in gene '+genes[i]+'\n')
							if num_snps > 0:
								winexclcount, file_count, win_count, num_snps = self.metric_calculation_gene(infile, outfile, winexclcount, file_count, win_count, num_snps)
								outkeep.write(genes[i]+'\n')
							else:
								print('No SNPs in candidate gene '+genes[i])
								noSNP_count +=1
								outfile2.write(genes[i]+'\n')
						outkeep.close()
						outfile2.close()
						logsnp.close()
				else:
					print('No genes in file.')
				if noSNP_count == 0:
					os.remove(outputdir+contrast+'_'+baseparname+'_NO_SNPS.txt')
			print('\nAnalyzed '+str(win_count)+' genes')#'+str(file_count)+' files')
			print('\t!! '+str(noSNP_count)+' genes had 0 SNPs!!\n')
			logfile.write('Analyzed '+str(win_count)+' genes.\n')
			logfile.write(str(noSNP_count)+' genes had 0 SNPs.\n')
		outfile.close()
		# logfile.close()


		# inherit objects to next step
		self.noSNP_count = noSNP_count
		self.baseparname = baseparname
		self.logfile = logfile


	def step4_gene(self):
		###### STEP 4 GENE ######
		## Remove genes with no SNPs from the gene interval bed/txt file and add DD, Tajimas D and Fai and Wu's H
		args = self.args
		outputdir = self.outputdir
		contrast = self.contrast
		noSNP_count = self.noSNP_count
		baseparname = self.baseparname
		logdir = self.logdir
		logfile = self.logfile

		print('\n\tSTEP 4 GENE: For genes with SNPs add DD, Pi per population, Tajimas D per pop, thetas\n')
		logfile.write('\n\tSTEP 4 GENE: For genes with SNPs add DD, Pi per population, Tajimas D per pop, thetas\n')

		if noSNP_count > 0:

		# read all lines of gene interval file
			with open(outputdir+contrast+'_'+baseparname+'_gene_interval.txt', 'r') as f:
				lines = f.readlines()
			f.close()
			# write lines except the ones that are in the no-snps file
			keep = open(outputdir+contrast+'_'+baseparname+'_gene_interval_keep.txt', 'w')
			with open(outputdir+contrast+'_'+baseparname+'_keep_genes.txt', 'r') as delete_line:
				for line in delete_line:
					delline = line.replace('\n', '')
					for line in lines:
						linekeep = line.split()
						if delline == linekeep[0]:
							keep.write(line)
			keep.close()

			# combine (paste) the gene interval file with the gene metrics file into one table
			pastecmd = open(outputdir+'paste.unix', 'w')
			pastecmd.write('paste '+outputdir+contrast+'_'+baseparname+'_gene_interval_keep.txt '+outputdir+contrast+'_temp_gene_metrics_'+baseparname+args.suf+'.txt'+' > '+outputdir+contrast+'_gene_metrics_'+baseparname+args.suf+'.txt')
			pastecmd.close()
			# Execute in unix
			cmd = open(outputdir+'paste.unix', 'r')
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid,0)[1]
			os.remove(outputdir+'paste.unix')

		else:
			# combine (paste) the gene interval file with the gene metrics file into one table
			pastecmd = open(outputdir+'paste.unix', 'w')
			pastecmd.write('paste '+outputdir+contrast+'_'+baseparname+'_gene_interval.txt '+outputdir+contrast+'_temp_gene_metrics_'+baseparname+args.suf+'.txt'+' > '+outputdir+contrast+'_gene_metrics_'+baseparname+args.suf+'.txt')
			pastecmd.close()
			# Execute in unix
			cmd = open(outputdir+'paste.unix', 'r')
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid,0)[1]
			os.remove(outputdir+'paste.unix')

		# Add DD residuals
		inlist = []
		for dirName, subdirList, fileList in os.walk(outputdir):
			for file in fileList:
				if file.startswith(contrast+'_gene_metrics_') and file.endswith(args.suf+'.txt'):
					inlist.append(dirName+'/'+file)
		inlist_sort = natsorted(inlist)

		print('Add DDresiduals')

		for file in inlist_sort:
			count = 0

			snp_file = pd.read_table(file, header=0)
			slope, intercept, r_value, p_value, std_err = stats.linregress(snp_file['absdiff'], snp_file[args.coh1+'_pi'])
			snp_file['prediction'] = intercept + (slope*snp_file['absdiff'])
			print('Uncorrected slope '+str(slope)+' uncorrected intercept '+str(intercept))
			logfile.write('For DD residuals:\nUncorrected slope '+str(slope)+' uncorrected intercept '+str(intercept)+'\n')

		# calculate slope correction using var(u)
			varu = snp_file.varu
			mean_absdiff = statistics.mean(snp_file.absdiff)
			mean_pi1 = statistics.mean(snp_file[args.coh1+'_pi'])
			var_absdiff = statistics.variance(snp_file.absdiff)
			cor_slope = statistics.mean((var_absdiff/(var_absdiff - varu)))*slope
			cor_intercept = mean_pi1-(mean_absdiff*cor_slope)

			print('\nVar(u) of '+str(file)+' is '+str(statistics.mean(varu)))
			print('Corrected slope '+str(cor_slope)+' corrected intercept '+str(cor_intercept)+'\n')
			logfile.write('Corrected slope '+str(cor_slope)+' corrected intercept '+str(cor_intercept)+'\n')
			logfile.write('Var(u)='+str(statistics.mean(varu))+'\n')

			snp_file['cor_prediction'] = cor_intercept + (cor_slope*snp_file.absdiff)
			snp_file['DD'] = snp_file[args.coh1+'_pi'] - snp_file.cor_prediction

		# calculate Tajima's D
			delta1 = snp_file[args.coh1+'_pi']-snp_file[args.coh1+'_thetaw']
			delta2 = snp_file[args.coh2+'_pi']-snp_file[args.coh2+'_thetaw']
			snp_file['tajimas_D_'+args.coh1] = delta1/statistics.stdev(delta1)
			snp_file['tajimas_D_'+args.coh2] = delta2/statistics.stdev(delta2)

		# calculate Fai and Wu's H
			delta1 = snp_file[args.coh1+'_pi']-snp_file[args.coh1+'_thetah']
			delta2 = snp_file[args.coh2+'_pi']-snp_file[args.coh2+'_thetah']
			snp_file['faiwus_H_'+args.coh1] = delta1/statistics.stdev(delta1)
			snp_file['faiwus_H_'+args.coh2] = delta2/statistics.stdev(delta2)

			snp_file.to_csv(outputdir+contrast+'_WG_metrics_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_'+str(int(args.ovlp))+'ol'+args.suf+'.txt',sep="\t", index=False)
			   
			count += 1
		   
		print('Processed '+str(count)+' files for DDresiduals')
		print('Processed '+str(count)+' files for Tajimas D')
		print('Processed '+str(count)+' files for Fai and Wus H\n')
		# Remove files that are no longer necessary
		os.remove(outputdir+contrast+'_temp_gene_metrics_'+baseparname+args.suf+'.txt')
		os.remove(outputdir+contrast+'_gene_metrics_'+baseparname+args.suf+'.txt')

		# inherit the following file
		self.snp_file = snp_file
		self.logfile = logfile

	def step5_gene(self):
		###### STEP 5 GENE ######
		# Combine snpEff summary table with per gene metrics table
		args = self.args
		outputdir = self.outputdir
		contrast = self.contrast
		contrast = args.coh1 + args.coh2
		snp_file = self.snp_file
		genes_bed_sorted = self.genes_bed_sorted
		noSNP_count = self.noSNP_count
		baseparname = self.baseparname
		logdir = self.logdir
		logfile = self.logfile

		print('\n\tSTEP 5 GENE: Annotate outlier SNPs using snpEff\n')
		logfile.write('\n\tSTEP 5 GENE: Annotate outlier SNPs using snpEff\n')

		# SNPeff annotation of target genes
		# test whether annotation has been done already
		if os.path.exists(outputdir+args.coh1+'_'+baseparname+'_ann'+args.suf+'.vcf') == False:
			# search for WG vcfs
			vcf_list = []
			for dirName, subdirList, fileList in os.walk(args.o):
				for file in fileList:
					if file.endswith('WG_'+args.coh1+'.vcf'):
						vcf_list.append(dirName+file)
			for vcf in vcf_list:
				print('Annotating vcf for ')
				
				keep_bed = open(outputdir+contrast+'_kept_genes.bed', 'w')
				# create bedfile for annotation intervalls
				with open(outputdir+contrast+'_'+baseparname+'_gene_interval_keep.txt', 'r') as tkeep:
					tkeep.readline()
					for line in tkeep:
						gene, scaffold, start, end, gene_length = line.split()
						keep_bed.write(scaffold+'\t'+start+'\t'+end+'\n')
				keep_bed.close()
				# annotate the candidate genes using SNPeff
				snpeffcmd = open(outputdir+'snpeff_targets.unix', 'w')
				snpeffcmd.write('snpeff -fi '+outputdir+contrast+'_kept_genes.bed -ud 0 -onlyProtein LyV2 '+vcf)
				snpeffcmd.write(' > '+outputdir+args.coh1+'_'+baseparname+'_ann'+args.suf+'.vcf')
				snpeffcmd.close()
				# execute in unix
				cmd = (open(outputdir+'snpeff_targets.unix', 'r'))
				p = subprocess.Popen(cmd, shell=True)
				sts = os.waitpid(p.pid, 0)[1]

				snpmv = open(outputdir+'snpmv.unix', 'w')
				snpmv.write('grep -v following snpEff_genes.txt > '+outputdir+args.coh1+'_snpEff_genes.txt')
				snpmv.close()
				# execute in unix
				cmd = (open(outputdir+'snpmv.unix', 'r'))
				p = subprocess.Popen(cmd, shell=True)
				sts = os.waitpid(p.pid, 0)[1]

				os.remove(outputdir+'snpeff_targets.unix')
				os.remove(outputdir+'snpmv.unix')

				print(str(vcf)+' annotated for !!\n')
				logfile.write(str(vcf)+' annotated using SNPeff.\n')
				# obtain snpEff_genes.txt file
				ann_file = pd.read_table(outputdir+args.coh1+'_snpEff_genes.txt', header=0)
				# print(ann_file.columns)
				# inner join of both files
				merged_inner = pd.merge(left=snp_file, right=ann_file, left_on='gene', right_on='GeneId')
				merged_inner.to_csv(outputdir+contrast+'_'+args.coh1+'_annotated_candidates'+args.suf+'.txt', sep="\t", index=False)

			logfile.close()
			
		else:
			# obtain snpEff_genes.txt file
			ann_file = pd.read_table(outputdir+args.coh1+'_snpEff_genes.txt', header=0)
			# print(ann_file.columns)
			# inner join of both files
			merged_inner = pd.merge(left=snp_file, right=ann_file, left_on='gene', right_on='GeneId')
			merged_inner.to_csv(outputdir+contrast+'_'+args.coh1+'_annotated_candidates'+args.suf+'.txt', sep="\t", index=False)


	# define the sequences of steps that should be run
	def step1to5(self):
		self.step1_gene()
		self.step2_gene()
		self.step3_gene()
		self.step4_gene()
		self.step5_gene()

	def step1to2(self):
		self.step1_gene()
		self.step2_gene()
			
if __name__ == '__main__':
	import os, sys
