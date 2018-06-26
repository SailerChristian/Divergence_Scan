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
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages


## create Class 
class G1():
    def __init__(self, AFs_table=None, three_metrics=None, args=None, natsorted=natsorted):
        self.AFs_table, self.three_metrics = None, None
        self.args = args
        self.natsorted = natsorted

    def step1(self):
        ###### STEP 1 ######
        ## Prepare input data for analysis
        # test where to place the outputdirectory
        args = self.args
        cwd = os.getcwd()
        if args.o is 'na':
            outputdir = str(args.coh1+args.coh2)
        else:
            outputdir = str(args.o+args.coh1+args.coh2)
            if os.path.exists(args.o) == False:
                os.mkdir(args.o)
        print('\nSetup ...')

        # check if output directory exists, create it if necessary
        if os.path.exists(outputdir) == False:
            os.mkdir(outputdir)
            print('\n\tCreated directory '+outputdir)
        contrast = args.coh1+args.coh2
        
        print('\n\tSearching input directory for *_raw.table')
        tincoh1 = []
        tincoh2 = []
        for dirName, subdirList, fileList in os.walk(args.i, topdown = False):
            for fname in fileList:
                if fname.endswith('_raw.table') and args.coh1 in fname:
                    tincoh1.append(dirName+fname)
                elif fname.endswith('_raw.table') and args.coh2 in fname:
                    tincoh2.append(dirName+fname)
        incoh1 = natsorted(tincoh1)
        incoh2 = natsorted(tincoh2)

        ## test to read header and determine the number of columns for pasting
        with open(incoh1[0], 'r') as infile:
            header = infile.readline()
            head = header.split()
            ncol = len(head)
            print('\t'+str(ncol)+' columns in input tables')
        infile.close()
                                    
        # paste cohort 1 & cohort 2 into a table next to each other & remove CHROM POS from second cohort in joint table
        # yields: CHROM POS AC AN AC AN
        print('\nSTEP 1: Paste contrast AC tables\n')
        for i in range(len(incoh1)):
            print('\tProcessing '+incoh1[i])
            pastecmd = open('paste_'+args.coh1+'_'+args.coh2+'.unix', 'w')
            pastecmd.write('paste ')
            pastecmd.write(incoh1[i]+' '+incoh2[i]+' | cut -f -'+str(ncol))
            for j in range(3,ncol+1):
                pastecmd.write(','+str(ncol*2-(ncol-j)))
            pastecmd.write(' > '+outputdir+'/'+contrast+'_scaf_'+str(i+1)+'_temp.table')
            pastecmd.close()

            # run in unix
            cmd = open('paste_'+args.coh1+'_'+args.coh2+'.unix', 'r')
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

        os.remove('paste_'+args.coh1+'_'+args.coh2+'.unix') # removes unix script file from folder

        # variables to be inherited
        self.outputdir = outputdir
        self.contrast = contrast


    def step2(self):
        ###### STEP 2 ######
        ## This section is Jeff's FixedDerivedAlleleCheck.py
        # search directory for output of above step to make new input list
        args = self.args
        outputdir = self.outputdir
        contrast = self.contrast
        tintable = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for table in fileList:
                if table.startswith(contrast) and table.endswith('_temp.table'):
                    tintable.append(dirName+'/'+table)
        intable = natsorted(tintable)
        
        print('\nSTEP 2: Remove sites that are fixed in both cohorts and are missing in one but not the other cohort\n')
        print('\tFound '+str(len(intable))+' input tables')

        # obtain ANmax for each cohort
        print('\n\tGathering allele count info:\n')
        if os.path.exists(outputdir+'/process_log/') == False:
            os.mkdir(outputdir+'/process_log/')
        logdir = outputdir+'/process_log/'
        self.logdir = logdir

        for i in range(len(intable)):
            with open(intable[i], 'r') as infile:
                infile.readline()
                tANmax_1 = []
                tANmax_2 = []
                for line in infile:
                    scaffold, position, AC_1, AN_1, AC_2, AN_2 = line.split()
                    tANmax_1.append(int(AN_1))
                    tANmax_2.append(int(AN_2))
                ANmax_1 = max(tANmax_1)
                ANmax_2 = max(tANmax_2)
                print('\tFor '+intable[i]+': ')
                print('\tANmax = '+str(ANmax_1)+' in cohort '+args.coh1)
                print('\tANmax = '+str(ANmax_2)+' in cohort '+args.coh2)

                # filter
                infile.seek(0)
                with open(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'_AC'+args.suf+'.table','w') as outfile:
                    header = infile.readline() # to read the first line of the file and thereby remove it from the rest of the loop
                    outfile.write('CHROM\tPOS\tAC1\tAN1\tAC2\tAN2\tvariant\n')

                    count = 0
                    target = 100000
                    Zero_count = 0
                    All_count = 0
                    non_variant_count = 0

                    for line in infile:
                        scaffold, position, AC_1, AN_1, AC_2, AN_2 = line.split()
                        if int(AN_1)>=int(round((1-args.per)*ANmax_1,0)) and int(AN_2)>=int(round((1-args.per)*ANmax_2,0)):
                            total_allele_count = int(AC_1)+int(AC_2)
                            total_alleles_sampled = int(AN_1)+int(AN_2)
                            # Count types of variants and add whether they are non-variant (NV) or SNPs (V)
                            if total_allele_count == total_alleles_sampled:
                                All_count +=1
                                outfile.write(scaffold+'\t'+position+'\t'+AC_1+'\t'+AN_1+'\t'+AC_2+'\t'+AN_2+'\t'+'NV\n')
                            elif total_allele_count == 0:
                                Zero_count +=1
                                outfile.write(scaffold+'\t'+position+'\t'+AC_1+'\t'+AN_1+'\t'+AC_2+'\t'+AN_2+'\t'+'NV\n')
                            elif int(AC_1) == int(AC_2):
                                non_variant_count += 1
                                outfile.write(scaffold+'\t'+position+'\t'+AC_1+'\t'+AN_1+'\t'+AC_2+'\t'+AN_2+'\t'+'NV\n')
                            else:
                                outfile.write(scaffold+'\t'+position+'\t'+AC_1+'\t'+AN_1+'\t'+AC_2+'\t'+AN_2+'\t'+'V\n')
                            count += 1

                    print('\tTotal sites passing filter in file: '+str(count))
                    print('\tSites with total allele count = 0: '+str(Zero_count))
                    print('\tVariants fixed in both cohorts: '+str(All_count))
                    print('\tNon-variant sites: '+str(non_variant_count))
                    print('\tVariant sites: '+str(count-Zero_count-All_count-non_variant_count)+'\n')
                    os.remove(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'_temp.table') # removes output files of part 1, as they are no longer needed

                    # write log file with counts
                    with open(logdir+contrast+'_scaf_'+str(i+1)+'_inputfilter.csv', 'w') as logfile:
                        logfile.write('ANmax in cohort '+args.coh1+','+str(ANmax_1)+'\n')
                        logfile.write('ANmax in cohort '+args.coh2+','+str(ANmax_2)+'\n')
                        logfile.write('Total sites passing in file,'+str(count)+'\n')
                        logfile.write('Sites with total allele count = 0,\t'+str(Zero_count)+'\n')
                        logfile.write('Variants fixed in both cohorts,'+str(All_count)+'\n')
                        logfile.write('Non-variant sites,'+str(non_variant_count)+'\n')
                        logfile.write('Variant sites (SNPs),'+str(count-Zero_count-All_count-non_variant_count))


    def step3(self):
        ###### STEP 3 ######
        ## jump highdepthFilter and remove DP column, Jeff's script adjusted
        args = self.args
        outputdir = self.outputdir
        contrast = self.contrast

        print('\nSTEP 3: Concatenate files into 1 WG file')
        # tintable = []
        # for dirName, subdirList, fileList in os.walk(outputdir):
        #     for table in fileList:
        #         if table.startswith(contrast) and table.endswith('.tab'):
        #             tintable.append(dirName+'/'+table)
        # intable = natsorted(tintable)

        # for i in range(len(intable)):
        #     with open(intable[i],'r') as infile:
        #         infile.readline()
        #         with open(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'_AFs'+args.suf+'.table','w') as outfile:
        #             outfile.write('CHROM\tPOS\tAC1\tAN1\tAC2\tAN2\n')

        #             for line in infile:
        #                 scaffold, position, AC_1, AN_1, DP_1, AC_2, AN_2, DP_2 = line.split()
        #                 outfile.write(scaffold+'\t'+position+'\t'+AC_1+'\t'+AN_1+'\t'+AC_2+'\t'+AN_2+'\n')
        #             os.remove(outputdir+'/'+contrast+'_scaf_'+str(i+1)+'.tab') # removes output files of part 2, as they are no longer needed
        # outfile.close()
        # infile.close()

        ## concatenate the AC tables to one WG (whole genome) table
        AFs = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith('_AC'+args.suf+'.table') == True:
                    AFs.append(dirName+'/'+file)
        AFs_sort = natsorted(AFs)

        with open(outputdir+'/'+contrast+'_WG_AC'+args.suf+'.txt','w') as outfile:
            for i in range(len(AFs_sort)):
                with open(AFs_sort[i], 'r') as infile:
                    infile.readline()
                    for line in infile:
                        outfile.write(line)
        outfile.close()
        infile.close()


    def metric_calculation(self, infile, outfile, winexclcount, file_count, win_count, num_sites):
        args = self.args
        scaf = []
        pos = []
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
        fst_4C = []
        dxy = []
        varu = []
        site_count = 0
        for i in range(num_sites):
            data = infile.readline()
            scaffold, position, AC_1, AN_1, AC_2, AN_2, variant = data.split()
            scaf.append(scaffold)
            pos.append(int(position))
            if int(AC_1) > 0 and int(AC_1) != int(AN_1):
                num_snps_1 += 1
            if int(AC_2) > 0 and int(AC_2) != int(AN_2):
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
            absdiff.append(abs(allele_c1_freq-allele_c2_freq))
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
            try:
                ht = 2*allele_freq*(1-allele_freq)
                h1 = 2*allele_c1_freq*(1-allele_c1_freq)
                h2 = 2*allele_c2_freq*(1-allele_c2_freq)
                h12 = ((h1*int(AN_1))+(h2*int(AN_2)))/(int(AN_1)+int(AN_2))
                allele_fst = abs(ht-h12)/ht
                fst.append(allele_fst)
            except ZeroDivisionError:
                pass
            # alternate heterozygosity based on autoteraploid frequencies 4p^3q+6p^2q^2+4pq^3
            try:
                ht_alt = 4*allele_freq**3*(1-allele_freq) + 6*allele_freq**2*(1-allele_freq)**2 + 4*allele_freq*(1-allele_freq)**3
                h1_alt = 4*allele_c1_freq**3*(1-allele_c1_freq) + 6*allele_c1_freq**2*(1-allele_c1_freq)**2 + 4*allele_c1_freq*(1-allele_c1_freq)**3
                h2_alt = 4*allele_c2_freq**3*(1-allele_c2_freq) + 6*allele_c2_freq**2*(1-allele_c2_freq)**2 + 4*allele_c2_freq*(1-allele_c2_freq)**3
                h12_alt = ((h1_alt*int(AN_1))+(h2_alt*int(AN_2)))/(int(AN_1)+int(AN_2))
                fst_alt = (ht_alt-h12_alt)/ht_alt
                fst_4C.append(fst_alt)
            except ZeroDivisionError:
                pass
            dxy.append((allele_c1_freq*(1-allele_c2_freq))+(allele_c2_freq*(1-allele_c1_freq)))
            varu1 = (1/(int(AN_1)*int(AN_1)))*allele_c1_freq*(1-allele_c1_freq)
            varu2 = (1/(int(AN_2)*int(AN_2)))*allele_c2_freq*(1-allele_c2_freq)
            varu.append(varu1 + varu2)

        wstart = min(pos)
        wend = max(pos)
        wlength = wend-wstart
        wmid = wstart+(wlength/2)
        mean_c1_freq = statistics.mean(c1_freq)
        mean_c2_freq = statistics.mean(c2_freq)
        mean_absdiff = statistics.mean(absdiff)
        mean_c1_minus_c2 = statistics.mean(c1_minus_c2)
        mean_c2_minus_c1 = statistics.mean(c2_minus_c1)
        a1 = sum_rec_1/args.snps
        a2 = sum_rec_2/args.snps
        sum_theta_pi1 = sum(theta_pi1)
        sum_theta_pi2 = sum(theta_pi2)
        pi1 = sum(theta_pi1)/len(pos)
        pi2 = sum(theta_pi2)/len(pos)
        theta_w1 = num_snps_1/a1
        theta_w2 = num_snps_2/a2
        sum_theta_h1 = sum(theta_h1)
        sum_theta_h2 = sum(theta_h2)
        mean_fst = statistics.mean(fst)
        mean_fst_4C = statistics.mean(fst_4C)
        sum_dxy = sum(dxy)
        win_dxy = (1/len(pos))*sum_dxy
        mean_varu = statistics.mean(varu)

        if wlength <= 26560:
            outfile.write(scaf[0]+'\t'+
                          str(wstart)+'\t'+
                          str(wend)+'\t'+
                          str(wmid)+'\t'+
                          str(wlength)+'\t'+
                          str(len(pos))+'\t'+
                          str(num_snps_1)+'\t'+
                          str(num_snps_2)+'\t'+
                          str(args.snps)+'\t'+
                          str(mean_c1_freq)+'\t'+
                          str(mean_c2_freq)+'\t'+
                          str(mean_absdiff)+'\t'+
                          str(mean_c1_minus_c2)+'\t'+
                          str(mean_c2_minus_c1)+'\t'+
                          str(mean_varu)+'\t'+
                          str(pi1)+'\t'+
                          str(pi2)+'\t'+
                          str(sum_theta_pi1)+'\t'+
                          str(sum_theta_pi2)+'\t'+
                          str(theta_w1)+'\t'+
                          str(theta_w2)+'\t'+
                          str(sum_theta_h1)+'\t'+
                          str(sum_theta_h2)+'\t'+
                          str(mean_fst)+'\t'+
                          str(mean_fst_4C)+'\t'+
                          str(win_dxy)+'\n')
        else:
            winexclcount +=1
        return winexclcount, file_count, win_count, num_sites


    def step4(self):
        ###### STEP 4 ######
        # calculate Dxy, Fst, pi, allele frequency and allele frequncy difference. Input required as:
        # scaffold, position, cohort1_allele_count, cohort1_all_alleles_count, cohort2_allele_count, cohort2_all_alleles_count
        args = self.args
        if args.o is 'na':
            outputdir = str(args.coh1+args.coh2)
        else:
            outputdir = str(args.o+args.coh1+args.coh2)

        contrast = args.coh1 + args.coh2
        self.outputdir = outputdir
        self.contrast = contrast
        print('\nSTEP 4: Calculate Dxy, Fst, DD')
        header = 'scaffold\tstart\tend\tmidpoint\tlength_bp\tnum_sites\tS_'+args.coh1+'\tS_'+args.coh2+'\tvariant_sites\t'+args.coh1+'_freq\t'+args.coh2+'_freq\tabsdiff\t'+args.coh1+'_'+args.coh2+'\t'+args.coh2+'_'+args.coh1+'\tvaru\t'+args.coh1+'_pi\t'+args.coh2+'_pi\t'+args.coh1+'_thetaPi\t'+args.coh2+'_thetaPi\t'+args.coh1+'_thetaw\t'+args.coh2+'_thetaw\t'+args.coh1+'_thetah\t'+args.coh2+'_thetah\tFst\tFst_4C\tDxy\n'
        with open(outputdir+'/'+contrast+'_metrics_WG_'+str(args.snps)+'SNPs'+args.suf+'.txt','w') as outfile:
            outfile.write(header)

            AFs = []
            for dirName, subdirList, fileList in os.walk(outputdir):
                for file in fileList:
                    if file.startswith(contrast) and file.endswith('_AC'+args.suf+'.table') == True:
                        AFs.append(dirName+'/'+file)
            AFs_sort = natsorted(AFs)

            print('\n\tFound '+str(len(AFs_sort))+' files starting with '+contrast+' and ending with AC'+args.suf+'.table\n')

            file_count = 0
            win_count = 0
            winexclcount = 0

            for AF in AFs_sort:
                count = 0 
                var_count = 0
                with open(AF,'r') as infile:
                    infile.readline() # to remove header from count
                    for line in infile:
                        count += 1 # counts the total number sites in file
                        # I cannot inherit the number from step 2 when changing args.snps. Redo the variant counting to get args.snps == num_var_sites per window
                        scaffold, position, AC_1, AN_1, AC_2, AN_2, variant = line.split()
                        total_allele_count = int(AC_1)+int(AC_2)
                        total_alleles_sampled = int(AN_1)+int(AN_2)
                        if total_allele_count == total_alleles_sampled:
                            pass
                        elif total_allele_count == 0:
                            pass
                        elif int(AC_1) == int(AC_2):
                            pass
                        else:
                            var_count += 1 # counts only SNPs
                    print('\tFound '+str(count)+' sites in file')
                    print('\tFound '+str(var_count)+' variant sites (SNPs) in file')
                    num_win_var = int(var_count/args.snps)
                    print('\tProcessing '+str(num_win_var)+' windows with '+str(args.snps)+' variants in '+AF)

                    # num_win = int(count/args.snps)
                    # print('\tProcessing '+str(num_win)+' windows in '+AF)

                    if num_win_var > 0:
                        file_count += 1
                        win_count += num_win_var

                        # with open(AF,'r') as infile:
                        #     infile.readline() # to read the header and start reading the numbers for calcluation from the second line of the file
                        #     for i in range(num_win_var):
                        #         winexclcount, file_count, win_count = self.metric_calculation(infile, outfile, winexclcount, file_count, win_count)

                    # Create a list with window lenghts that each contain args.snps variant sites.
                        with(open(AF, 'r')) as infile:
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
                                    scaffold, position, AC_1, AN_1, AC_2, AN_2, variant = line.split()
                                    if variant == 'V':
                                        snp_count += 1
                                    if snp_count%target == 0:
                                        # print(str(snp_count)+' SNPs '+str(site_count)+' sites')
                                        break
                                win_len.append(site_count)

                            else:
                                print(len(win_len))
                            infile.seek(0) # go back to the start
                            infile.readline() # read header
                            for i in range(num_win_var):
                                num_sites = int(win_len[i])
                                winexclcount, file_count, win_count, num_sites = self.metric_calculation(infile, outfile, winexclcount, file_count, win_count, num_sites)

                            infile.close()

        print('\nAnalyzed '+str(win_count)+' windows in '+str(file_count)+' files')
        print('Excluded '+str(winexclcount)+' windows longer than 26560bp\n')

        ###### STEP 4b ######
        # Add DD residuals
        inlist = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs'+args.suf+'.txt'):
                    inlist.append(dirName+'/'+file)
        inlist_sort = natsorted(inlist)

        print('\tAdd DDresiduals')

        for file in inlist_sort:
            count = 0

            snp_file = pd.read_table(file, header=0)
            slope, intercept, r_value, p_value, std_err = stats.linregress(snp_file['absdiff'], snp_file[args.coh1+'_pi'])
            snp_file['prediction'] = intercept + (slope*snp_file['absdiff'])
            print('Uncorrected slope '+str(slope)+' uncorrected intercept '+str(intercept))

        # calculate slope correction using var(u)
            varu = snp_file.varu
            mean_absdiff = statistics.mean(snp_file.absdiff)
            mean_pi1 = statistics.mean(snp_file[args.coh1+'_pi'])
            var_absdiff = statistics.variance(snp_file.absdiff)
            cor_slope = statistics.mean((var_absdiff/(var_absdiff - varu)))*slope
            cor_intercept = mean_pi1-(mean_absdiff*cor_slope)

            print('\nVar(u) of '+str(file)+' is '+str(statistics.mean(varu)))
            print('Corrected slope '+str(cor_slope)+' corrected intercept '+str(cor_intercept))
            
            snp_file['cor_prediction'] = cor_intercept + (cor_slope*snp_file.absdiff)
            snp_file['DD'] = snp_file[args.coh1+'_pi'] - snp_file.cor_prediction

        # calculate Tajima's D
            delta1 = snp_file[args.coh1+'_thetaPi']-snp_file[args.coh1+'_thetaw']
            delta2 = snp_file[args.coh2+'_thetaPi']-snp_file[args.coh2+'_thetaw']
            snp_file['tajimas_D_'+args.coh1] = delta1/statistics.stdev(delta1)
            snp_file['tajimas_D_'+args.coh2] = delta2/statistics.stdev(delta2)

            snp_file.to_csv(outputdir+'/'+contrast+'_WG_'+str(args.snps)+'SNPs_3metrics'+args.suf+'.txt',sep="\t", index=False)
            
            # count += 1
            
        # calculate Fai and Wu's H
            delta1 = snp_file[args.coh1+'_thetaPi']-snp_file[args.coh1+'_thetah']
            delta2 = snp_file[args.coh2+'_thetaPi']-snp_file[args.coh2+'_thetah']
            snp_file['faiwus_H_'+args.coh1] = delta1/statistics.stdev(delta1)
            snp_file['faiwus_H_'+args.coh2] = delta2/statistics.stdev(delta2)

            snp_file.to_csv(outputdir+'/'+contrast+'_WG_'+str(args.snps)+'SNPs_3metrics'+args.suf+'.txt',sep="\t", index=False)
            
            count += 1
        
        print('Processed '+str(count)+' files for DDresiduals')
        print('Processed '+str(count)+' files for Tajimas D')
        print('Processed '+str(count)+' files for Fai and Wus H')
        os.remove(outputdir+'/'+contrast+'_metrics_WG_'+str(args.snps)+'SNPs'+args.suf+'.txt')


    def step5(self):
        ###### STEP 5 ######
        # define outliers using top percentile cut-off
        print('\n\tSTEP 5: Select outlier windows')
        args = self.args
        if args.o is 'na':
            outputdir = str(args.coh1+args.coh2)
        else:
            outputdir = str(args.o+args.coh1+args.coh2)
        logdir = outputdir+'/process_log/'
        cwd = os.getcwd()
        if args.o is 'na':
            outputdir = str(args.coh1+args.coh2)
        else:
            outputdir = str(args.o+args.coh1+args.coh2)

        contrast = args.coh1 + args.coh2
        self.outputdir = outputdir
        self.contrast = contrast
        self.cwd = cwd

        out = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_3metrics'+args.suf+'.txt') == True:
                    out.append(dirName+'/'+file)
        out_sort = natsorted(out)

        for out in out_sort:
            with open(out, 'r') as infile:
                header = infile.readline()
                absafd = []
                afd_c1_c2 = []
                afd_c2_c1 = []
                pi_c1 = []
                pi_c2 = []
                thetaw_c1 = []
                thetaw_c2 = []
                thetah_c1 = []
                thetah_c2 = []
                dxy = []
                fst = []
                fst_4C = []
                dd = []
                tajD_c1 = []
                tajD_c2 = []
                tajH_c1 = []
                tajH_c2 = []
                faiwuH_c1 = []
                faiwuH_c2 = []
                wlength = []
                
                for line in infile:
                    scaffold, start, end, midpoint, length, num_sites, S_coh1, S_coh2, num_SNPs, coh1_freq, coh2_freq, absdiff, coh1_coh2, coh2_coh1, varu, coh1_pi, coh2_pi, coh1_thetaPi, coh2_thetaPi, coh1_thetaw, coh2_thetaw, coh1_thetah, coh2_thetah, Fst, Fst_4C, Dxy, prediction, cor_prediction, DD, tajimas_D_coh1, tajimas_D_coh2, faiwus_H_coh1, faiwus_H_coh2 = line.split()
                    absafd.append(float(absdiff))
                    afd_c1_c2.append(float(coh1_coh2))
                    afd_c2_c1.append(float(coh2_coh1))
                    pi_c1.append(float(coh1_pi))
                    pi_c2.append(float(coh2_pi))
                    thetaw_c1.append(float(coh1_thetaw))
                    thetaw_c2.append(float(coh2_thetaw))
                    thetah_c1.append(float(coh1_thetah))
                    thetah_c2.append(float(coh2_thetah))
                    dxy.append(float(Dxy))
                    fst.append(float(Fst))
                    fst_4C.append(float(Fst_4C))
                    dd.append(float(DD))
                    wlength.append(int(length))
                    tajD_c1.append(float(tajimas_D_coh1))
                    tajD_c2.append(float(tajimas_D_coh2))
                    faiwuH_c1.append(float(faiwus_H_coh1))
                    faiwuH_c2.append(float(faiwus_H_coh2))

                mean_afd = statistics.mean(afd_c1_c2)
                median_afd = statistics.mean(afd_c1_c2)
                sd_afd = statistics.stdev(afd_c1_c2)
                CV_afd = sd_afd/mean_afd
                out_afd_up = np.percentile(afd_c1_c2, (100-args.cut))
                out_afd_low = np.percentile(afd_c1_c2, args.cut)
                mean_abs = statistics.mean(absafd)
                median_abs = statistics.mean(absafd)
                sd_abs = statistics.stdev(absafd)
                CV_abs = sd_abs/mean_abs
                out_abs_up = np.percentile(absafd, (100-args.cut))
                out_abs_low = np.percentile(absafd, args.cut)
                mean_Dxy = statistics.mean(dxy)
                median_Dxy = statistics.median(dxy)
                sd_Dxy = statistics.stdev(dxy)
                CV_Dxy = sd_Dxy/mean_Dxy
                out_dxy = np.percentile(dxy, (100-args.cut))
                mean_Fst = statistics.mean(fst)
                median_Fst = statistics.median(fst)
                sd_Fst = statistics.stdev(fst)
                CV_Fst = sd_Fst/mean_Fst
                out_fst = np.percentile(fst, (100-args.cut))
                mean_Fst_4C = statistics.mean(fst_4C)
                median_Fst_4C = statistics.median(fst_4C)
                sd_Fst_4C = statistics.stdev(fst_4C)
                CV_Fst_4C = sd_Fst_4C/mean_Fst_4C
                out_fst_4C = np.percentile(fst_4C, (100-args.cut))
                mean_DD = statistics.mean(dd)
                median_DD = statistics.median(dd)
                sd_DD = statistics.stdev(dd)
                CV_DD = sd_DD/mean_DD
                out_dd = np.percentile(dd, args.cut)
                mean_win = statistics.mean(wlength)
                median_win = statistics.median(wlength)
                sd_win = statistics.stdev(wlength)
                CV_win = sd_win/mean_win
                out_win = np.percentile(wlength, (100-args.cut))
                mean_tajd1 = statistics.mean(tajD_c1)
                median_tajd1 = statistics.median(tajD_c1)
                sd_tajd1 = statistics.stdev(tajD_c1)
                CV_tajd1 = sd_tajd1/mean_tajd1
                out_tajd1 = np.percentile(tajD_c1, (100-args.cut))
                mean_tajd2 = statistics.mean(tajD_c2)
                median_tajd2 = statistics.median(tajD_c2)
                sd_tajd2 = statistics.stdev(tajD_c2)
                CV_tajd2 = sd_tajd2/mean_tajd2
                out_tajd2 = np.percentile(tajD_c2, (100-args.cut))
                mean_faiwu1 = statistics.mean(faiwuH_c1)
                median_faiwu1 = statistics.median(faiwuH_c1)
                sd_faiwu1 = statistics.stdev(faiwuH_c1)
                CV_faiwu1 = sd_faiwu1/mean_faiwu1
                out_faiwu1 = np.percentile(faiwuH_c1, (100-args.cut))
                mean_faiwu2 = statistics.mean(faiwuH_c2)
                median_faiwu2 = statistics.median(faiwuH_c2)
                sd_faiwu2 = statistics.stdev(faiwuH_c2)
                CV_faiwu2 = sd_faiwu2/mean_faiwu2
                out_faiwu2 = np.percentile(faiwuH_c2, (100-args.cut))

                outlier_values = open(outputdir+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_descriptive_stats'+args.suf+'.txt', 'w')
                outlier_values.write('stat\tAbsDiff\tAFD_1_2\tDxy\tFst\tFst_4C\tDD\twin_length\tTajD_'+args.coh1+'\tTajD_'+args.coh2+'\tFaiWusH_'+args.coh1+'\tFaiWusH_'+args.coh2+'\n'+
                                     'min\t'+str(round(min(absafd),4))+'\t'+str(round(min(afd_c1_c2),4))+'\t'+str(round(min(dxy),4))+'\t'+str(round(min(fst),4))+'\t'+str(round(min(fst_4C),4))+'\t'+str(round(min(dd),4))+'\t'+str(round(min(wlength),0))+'\t'+str(round(min(tajD_c1),4))+'\t'+str(round(min(tajD_c2),4))+'\t'+str(round(min(faiwuH_c1),4))+'\t'+str(round(min(faiwuH_c2),4))+'\n'+
                                     'mean\t'+str(round(mean_abs,4))+'\t'+str(round(mean_afd,4))+'\t'+str(round(mean_Dxy,4))+'\t'+str(round(mean_Fst,4))+'\t'+str(round(mean_Fst_4C,4))+'\t'+str(round(mean_DD,4))+'\t'+str(round(mean_win,0))+'\t'+str(round(mean_tajd1,4))+'\t'+str(round(mean_tajd2,0))+'\t'+str(round(mean_faiwu1,4))+'\t'+str(round(mean_faiwu2,0))+'\n'+
                                     'median\t'+str(round(median_abs,4))+'\t'+str(round(median_afd,4))+'\t'+str(round(median_Dxy,4))+'\t'+str(round(median_Fst,4))+'\t'+str(round(median_Fst_4C,4))+'\t'+str(round(median_DD,4))+'\t'+str(round(median_win,0))+'\t'+str(round(median_tajd1,4))+'\t'+str(round(median_tajd2,0))+'\t'+str(round(median_faiwu1,4))+'\t'+str(round(median_faiwu2,0))+'\n'+
                                     'sd\t'+str(round(sd_abs,4))+'\t'+str(round(sd_afd,4))+'\t'+str(round(sd_Dxy,4))+'\t'+str(round(sd_Fst,4))+'\t'+str(round(sd_Fst_4C,4))+'\t'+str(round(sd_DD,4))+'\t'+str(round(sd_win,0))+'\t'+str(round(sd_tajd1,4))+'\t'+str(round(sd_tajd2,0))+'\t'+str(round(sd_faiwu1,4))+'\t'+str(round(sd_faiwu2,0))+'\n'+
                                     'CV\t'+str(round(CV_abs,4))+'\t'+str(round(CV_afd,4))+'\t'+str(round(CV_Dxy,4))+'\t'+str(round(CV_Fst,4))+'\t'+str(round(CV_Fst_4C,4))+'\t'+str(round(CV_DD,4))+'\t'+str(round(CV_win,0))+'\t'+str(round(CV_tajd1,4))+'\t'+str(round(CV_tajd2,0))+'\t'+str(round(CV_faiwu1,4))+'\t'+str(round(CV_faiwu2,0))+'\n'+
                                     'max\t'+str(round(max(absafd),4))+'\t'+str(round(max(afd_c1_c2),4))+'\t'+str(round(max(dxy),4))+'\t'+str(round(max(fst),4))+'\t'+str(round(max(fst_4C),4))+'\t'+str(round(max(dd),4))+'\t'+str(round(max(wlength),0))+'\t'+str(round(max(tajD_c1),4))+'\t'+str(round(max(tajD_c2),0))+'\t'+str(round(max(faiwuH_c1),4))+'\t'+str(round(max(faiwuH_c2),0))+'\n'+
                                     'Cutoff\t'+str(round(out_abs_up,4))+'\t'+str(round(out_afd_up,4))+'\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_fst_4C,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\t'+str(round(out_tajd1,4))+'\t'+str(round(out_tajd2,0))+'\t'+str(round(out_faiwu1,4))+'\t'+str(round(out_faiwu2,0))+'\n'+
                                     'Cutoff2\t'+str(round(out_abs_low,4))+'\t'+str(round(out_afd_low,4))+'\t'+str(round(out_dxy,4))+'\t'+str(round(out_fst,4))+'\t'+str(round(out_fst_4C,4))+'\t'+str(round(out_dd,4))+'\t'+str(round(out_win,0))+'\t'+str(round(out_tajd1,4))+'\t'+str(round(out_tajd2,0))+'\t'+str(round(out_faiwu1,4))+'\t'+str(round(out_faiwu2,0))+'\n')


                # use pandas to read infiles as dataframe and select outlier combinations
                c1_c2 = str(args.coh1+'_'+args.coh2)
                df_dat = pd.read_table(out)
                df_dat['afdout'] = np.where(df_dat[c1_c2] >= out_afd_up, 1, 0)
                df_dat['afd21out'] = np.where(df_dat[c1_c2] <= out_afd_low, 1, 0)
                df_dat['dxyout'] = np.where(df_dat.Dxy >= out_dxy, 1, 0)
                if args.ploidy == 2:
                    df_dat['fstout'] = np.where(df_dat.Fst >= out_fst, 1, 0)
                else:
                    df_dat['fstout'] = np.where(df_dat.Fst_4C >= out_fst_4C, 1, 0)
                df_dat['ddout'] = np.where(df_dat.DD <= out_dd, 1, 0)
                # doubles
                df_dat['afddxy'] = df_dat.afdout+df_dat.dxyout
                df_dat['afdfst'] = df_dat.afdout+df_dat.fstout
                df_dat['afddd'] = df_dat.afdout+df_dat.ddout
                df_dat['afd21dxy'] = df_dat.afd21out+df_dat.dxyout
                df_dat['afd21fst'] = df_dat.afd21out+df_dat.fstout
                df_dat['afd21dd'] = df_dat.afd21out+df_dat.ddout
                df_dat['dxyfst'] = df_dat.dxyout+df_dat.fstout
                df_dat['dxydd'] = df_dat.dxyout+df_dat.ddout
                df_dat['fstdd'] = df_dat.fstout+df_dat.ddout
                # triples
                df_dat['afddxyfst'] = df_dat.afdout+df_dat.dxyout+df_dat.fstout
                df_dat['afddxydd'] = df_dat.afdout+df_dat.dxyout+df_dat.ddout
                df_dat['afdfstdd'] = df_dat.afdout+df_dat.fstout+df_dat.ddout
                df_dat['afd21dxyfst'] = df_dat.afd21out+df_dat.dxyout+df_dat.fstout
                df_dat['afd21dxydd'] = df_dat.afd21out+df_dat.dxyout+df_dat.ddout
                df_dat['afd21fstdd'] = df_dat.afd21out+df_dat.fstout+df_dat.ddout
                df_dat['dxyfstdd'] = df_dat.dxyout+df_dat.fstout+df_dat.ddout
                # quadruple
                df_dat['afddxyfstdd'] = df_dat.afdout+df_dat.dxyout+df_dat.fstout+df_dat.ddout
                df_dat['afd21dxyfstdd'] = df_dat.afd21out+df_dat.dxyout+df_dat.fstout+df_dat.ddout

                df_dat['bedstart'] = df_dat.start-1
                print('\nFor top '+str(args.cut)+' percent we find:\n')
                # print(str(sum(1 for x in df_dat.afdout if x == 1))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4)))
                # print(str(sum(1 for x in df_dat.afd21out if x == 1))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4)))
                print(str(sum(1 for x in df_dat.dxyout if x == 1))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4)))
                if args.ploidy == 2:
                    print(str(sum(1 for x in df_dat.fstout if x == 1))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4)))
                else:
                    print(str(sum(1 for x in df_dat.fstout if x == 1))+'\toutlier windows for Fst_4C values of >='+str(round(out_fst_4C, 4)))                   
                print(str(sum(1 for x in df_dat.ddout if x == 1))+'\toutlier windows for DD values of <='+str(round(out_dd, 4)))
                print('\tof which are double outliers for:')
                # print(str(sum(1 for x in df_dat.afddxy if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and Dxy values of >='+str(round(out_dxy, 4)))
                # if args.ploidy == 2:
                #     print(str(sum(1 for x in df_dat.afdfst if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and Fst values of >='+str(round(out_fst, 4)))
                # else:
                #     print(str(sum(1 for x in df_dat.afdfst if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and Fst_4C values of >='+str(round(out_fst_4C, 4)))
                # print(str(sum(1 for x in df_dat.afddd if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and DD values of >='+str(round(out_dd, 4)))
                # print(str(sum(1 for x in df_dat.afd21dxy if x == 2))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+' and Dxy values of >='+str(round(out_dxy, 4)))
                # if args.ploidy == 2:
                #     print(str(sum(1 for x in df_dat.afd21fst if x == 2))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+' and Fst values of >='+str(round(out_fst, 4)))
                # else:
                #     print(str(sum(1 for x in df_dat.afd21fst if x == 2))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+' and Fst_4C values of >='+str(round(out_fst_4C, 4)))
                # print(str(sum(1 for x in df_dat.afd21dd if x == 2))+'\toutlier windows for AFD_2_1 values of <='+str(round(out_afd_low, 4))+' and DD values of >='+str(round(out_dd, 4)))
                if args.ploidy == 2:
                    print(str(sum(1 for x in df_dat.dxyfst if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of >='+str(round(out_fst, 4)))
                else:
                    print(str(sum(1 for x in df_dat.dxyfst if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and Fst_4C values of >='+str(round(out_fst_4C, 4)))
                print(str(sum(1 for x in df_dat.dxydd if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
                if args.ploidy == 2:
                    print(str(sum(1 for x in df_dat.fstdd if x == 2))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
                else:
                    print(str(sum(1 for x in df_dat.fstdd if x == 2))+'\toutlier windows for Fst_4C values of >='+str(round(out_fst_4C, 4))+' and DD values of <='+str(round(out_dd, 4)))
                print('\tof which are triple outliers for:')
                # if args.ploidy == 2:
                #     print(str(sum(1 for x in df_dat.afddxyfst if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of <='+str(round(out_fst, 4)))
                # else:
                #     print(str(sum(1 for x in df_dat.afddxyfst if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst_4C values of <='+str(round(out_fst_4C, 4)))
                # print(str(sum(1 for x in df_dat.afddxydd if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
                # if args.ploidy == 2:
                #     print(str(sum(1 for x in df_dat.afdfstdd if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
                # else:
                #     print(str(sum(1 for x in df_dat.afdfstdd if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Fst_4C values of >='+str(round(out_fst_4C, 4))+' and DD values of <='+str(round(out_dd, 4)))
                # if args.ploidy == 2:
                #     print(str(sum(1 for x in df_dat.afd21dxyfst if x == 3))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of <='+str(round(out_fst, 4)))
                # else:
                #     print(str(sum(1 for x in df_dat.afd21dxyfst if x == 3))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst_4C values of <='+str(round(out_fst_4C, 4)))
                # print(str(sum(1 for x in df_dat.afd21dxydd if x == 3))+'\toutlier windows for AFD_2_1 values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4)))
                # print(str(sum(1 for x in df_dat.afd21fstdd if x == 3))+'\toutlier windows for AFD_2_1 values of <='+str(round(out_afd_low, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
                if args.ploidy == 2:
                    print(str(sum(1 for x in df_dat.dxyfstdd if x == 3))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4)))
                else:
                    print(str(sum(1 for x in df_dat.dxyfstdd if x == 3))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+', Fst_4C values of >='+str(round(out_fst_4C, 4))+' and DD values of <='+str(round(out_dd, 4)))
                print('\tof which are quadruple outliers for:')
                # if args.ploidy == 2:
                #     print(str(sum(1 for x in df_dat.afddxyfstdd if x == 4))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                # else:
                #     print(str(sum(1 for x in df_dat.afddxyfstdd if x == 4))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+', Fst_4C values of >='+str(round(out_fst_4C, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                # print(str(sum(1 for x in df_dat.afd21dxyfstdd if x == 4))+'\toutlier windows for AFD_2_1 values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')

                # write log file with counts
                with open(logdir+contrast+'_top_'+str(args.cut)+'_percent_outlier.txt', 'w') as logfile:
                    # logfile.write(str(sum(1 for x in df_dat.afdout if x == 1))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+'\n')
                    # logfile.write(str(sum(1 for x in df_dat.afd21out if x == 1))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+'\n')
                    logfile.write(str(sum(1 for x in df_dat.dxyout if x == 1))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+'\n')
                    if args.ploidy == 2:
                        logfile.write(str(sum(1 for x in df_dat.fstout if x == 1))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4))+'\n')
                    else:
                        logfile.write(str(sum(1 for x in df_dat.fstout if x == 1))+'\toutlier windows for Fst_4C values of >='+str(round(out_fst_4C, 4))+'\n')                   
                    logfile.write(str(sum(1 for x in df_dat.ddout if x == 1))+'\toutlier windows for DD values of <='+str(round(out_dd, 4))+'\n')
                    logfile.write('\tof which are double outliers for:\n')
                    # logfile.write(str(sum(1 for x in df_dat.afddxy if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and Dxy values of >='+str(round(out_dxy, 4))+'\n')
                    # if args.ploidy == 2:
                    #     logfile.write(str(sum(1 for x in df_dat.afdfst if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and Fst values of >='+str(round(out_fst, 4))+'\n')
                    # else:
                    #     logfile.write(str(sum(1 for x in df_dat.afdfst if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and Fst_4C values of >='+str(round(out_fst_4C, 4))+'\n')
                    # logfile.write(str(sum(1 for x in df_dat.afddd if x == 2))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+' and DD values of >='+str(round(out_dd, 4))+'\n')
                    # logfile.write(str(sum(1 for x in df_dat.afd21dxy if x == 2))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+' and Dxy values of >='+str(round(out_dxy, 4))+'\n')
                    if args.ploidy == 2:
                        # logfile.write(str(sum(1 for x in df_dat.afd21fst if x == 2))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+' and Fst values of >='+str(round(out_fst, 4))+'\n')
                        logfile.write(str(sum(1 for x in df_dat.dxyfst if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of >='+str(round(out_fst, 4))+'\n')
                        logfile.write(str(sum(1 for x in df_dat.fstdd if x == 2))+'\toutlier windows for Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                    else:
                        # logfile.write(str(sum(1 for x in df_dat.afd21fst if x == 2))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+' and Fst_4C values of >='+str(round(out_fst_4C, 4))+'\n')
                        logfile.write(str(sum(1 for x in df_dat.dxyfst if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and Fst_4C values of >='+str(round(out_fst_4C, 4))+'\n')
                        logfile.write(str(sum(1 for x in df_dat.fstdd if x == 2))+'\toutlier windows for Fst_4C values of >='+str(round(out_fst_4C, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                    logfile.write(str(sum(1 for x in df_dat.dxydd if x == 2))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                    logfile.write('\tof which are triple outliers for:\n')
                    if args.ploidy == 2:
                        # logfile.write(str(sum(1 for x in df_dat.afddxyfst if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of <='+str(round(out_fst, 4))+'\n')
                        # logfile.write(str(sum(1 for x in df_dat.afdfstdd if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                        # logfile.write(str(sum(1 for x in df_dat.afd21dxyfst if x == 3))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst values of <='+str(round(out_fst, 4))+'\n')
                        logfile.write(str(sum(1 for x in df_dat.dxyfstdd if x == 3))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                    else:
                        # logfile.write(str(sum(1 for x in df_dat.afddxyfst if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst_4C values of <='+str(round(out_fst_4C, 4))+'\n')
                        # logfile.write(str(sum(1 for x in df_dat.afdfstdd if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Fst_4C values of >='+str(round(out_fst_4C, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                        # logfile.write(str(sum(1 for x in df_dat.afd21dxyfst if x == 3))+'\toutlier windows for AFD('+args.coh2+'-'+args.coh1+') values of <='+str(round(out_afd_low, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and Fst_4C values of <='+str(round(out_fst_4C, 4))+'\n')
                        logfile.write(str(sum(1 for x in df_dat.dxyfstdd if x == 3))+'\toutlier windows for Dxy values of >='+str(round(out_dxy, 4))+', Fst_4C values of >='+str(round(out_fst_4C, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                    # logfile.write(str(sum(1 for x in df_dat.afddxydd if x == 3))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                    # logfile.write('\tof which are quadruple outliers for:\n')
                    # if args.ploidy == 2:
                    #     logfile.write(str(sum(1 for x in df_dat.afddxyfstdd if x == 4))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+', Fst values of >='+str(round(out_fst, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                    # else:
                    #     logfile.write(str(sum(1 for x in df_dat.afddxyfstdd if x == 4))+'\toutlier windows for AFD('+args.coh1+'-'+args.coh2+') values of >='+str(round(out_afd_up, 4))+', Dxy values of >='+str(round(out_dxy, 4))+', Fst_4C values of >='+str(round(out_fst_4C, 4))+' and DD values of <='+str(round(out_dd, 4))+'\n')
                logfile.close()

                print('Select oulier windows\n')

                file_basename = cwd+'/'+outputdir+'/'+args.coh1+args.coh2+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_'
                df_dat.to_csv(file_basename+'allsites'+args.suf+'.csv', index=False)
                # select all windows that are outliers for at least one metric
                df_outlier = df_dat[((df_dat.dxyout != 0) | (df_dat.fstout != 0) | (df_dat.ddout != 0))]
                df_outlier.to_csv(file_basename+'ALL_outliers'+args.suf+'.csv', index=False)
                # write bedfile for gene retrieval and orthologous gene onthology match
                header = ["scaffold", "bedstart", "end"]
                df_outlier.to_csv(file_basename+'ALL_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # create inclusive outlier window lists, that is windows that are outliers for 2 or more metrics are listed in the 2 or more list as well as in the single metric lists
                # select AFD_1_2 only outliers
                # df_afd12 = df_dat[(df_dat.afdout == 1)]
                # if df_afd12.empty:
                #     print('No single AFD('+args.coh1+'-'+args.coh2+') outlier windows')
                # else:
                #     df_afd12.to_csv(file_basename+'afd12_outliers'+args.suf+'.csv', index=False)
                #     # write bedfile for gene retrieval and orthologous gene onthology match
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afd12.to_csv(file_basename+'afd12_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select AFD_2_1 only outliers
                # df_afd21 = df_dat[(df_dat.afd21out == 1)]
                # if df_afd21.empty:
                #     print('No single AFD('+args.coh2+'-'+args.coh1+') outlier windows')
                # else:
                #     df_afd21.to_csv(file_basename+'afd21_outliers'+args.suf+'.csv', index=False)
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afd21.to_csv(file_basename+'afd21_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select Dxy only outliers
                df_dxy = df_dat[(df_dat.dxyout == 1)] # exclusive df_dat[((df_dat.dxyout == 1) & (df_dat.fstout == 0) & (df_dat.ddout == 0))]
                if df_dxy.empty:
                    print('No single Dxy outlier windows')
                else:
                    df_dxy.to_csv(file_basename+'Dxy_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dxy.to_csv(file_basename+'Dxy_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select Fst only outliers
                df_fst = df_dat[(df_dat.fstout == 1)]
                if df_fst.empty:
                    print('No single Fst outlier windows')
                else:
                    df_fst.to_csv(file_basename+'Fst_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_fst.to_csv(file_basename+'Fst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select DD only outliers
                df_dd = df_dat[(df_dat.ddout == 1)]
                if df_dd.empty:
                    print('No single DD outlier windows')
                else:
                    df_dd.to_csv(file_basename+'DD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dd.to_csv(file_basename+'DD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select AFD_1_2 and Dxy double outliers
                # df_afddxy = df_dat[(df_dat.afddxy == 2) ]
                # if df_afddxy.empty:
                #     print('No double AFD('+args.coh1+'-'+args.coh2+') Dxy outlier windows')
                # else:
                #     df_afddxy.to_csv(file_basename+'afd12Dxy_outliers'+args.suf+'.csv', index=False)
                #     # write bedfile for gene retrieval and orthologous gene onthology match
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afddxy.to_csv(file_basename+'afd12Dxy_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # # select AFD_1_2 and Fst double outliers
                # df_afdfst = df_dat[(df_dat.afdfst == 2) ]
                # if df_afdfst.empty:
                #     print('No double AFD('+args.coh1+'-'+args.coh2+') Fst outlier windows')
                # else:
                #     df_afdfst.to_csv(file_basename+'afd12Fst_outliers'+args.suf+'.csv', index=False)
                #     # write bedfile for gene retrieval and orthologous gene onthology match
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afdfst.to_csv(file_basename+'afd12Fst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # # select AFD_1_2 and DD double outliers
                # df_afddd = df_dat[(df_dat.afddd == 2) ]
                # if df_afddd.empty:
                #     print('No double AFD('+args.coh1+'-'+args.coh2+') outlier windows')
                # else:
                #     df_afddd.to_csv(file_basename+'afd12DD_outliers'+args.suf+'.csv', index=False)
                #     # write bedfile for gene retrieval and orthologous gene onthology match
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afddd.to_csv(file_basename+'afd12DD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # # select AFD_2_1 and Dxy double outliers
                # df_afd21dxy = df_dat[(df_dat.afd21dxy == 2) ]
                # if df_afd21dxy.empty:
                #     print('No double AFD('+args.coh2+'-'+args.coh1+') Dxy outlier windows')
                # else:
                #     df_afd21dxy.to_csv(file_basename+'afd21Dxy_outliers'+args.suf+'.csv', index=False)
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afd21dxy.to_csv(file_basename+'afd21Dxy_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # # select AFD_2_1 and Fst double outliers
                # df_afd21fst = df_dat[(df_dat.afd21fst == 2) ]
                # if df_afd21fst.empty:
                #     print('No double AFD('+args.coh2+'-'+args.coh1+') Fst outlier windows')
                # else:
                #     df_afd21fst.to_csv(file_basename+'afd21Fst_outliers'+args.suf+'.csv', index=False)
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afd21fst.to_csv(file_basename+'afd21Fst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # select Dxy and Fst double outliers
                df_dxyfst = df_dat[(df_dat.dxyfst == 2) ]
                if df_dxyfst.empty:
                    print('No double Dxy Fst outlier windows')
                else:
                    df_dxyfst.to_csv(file_basename+'DxyFst_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dxyfst.to_csv(file_basename+'DxyFst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select Dxy and DD double outliers
                df_dxydd = df_dat[(df_dat.dxydd == 2)]
                if df_dxydd.empty:
                    print('No double Dxy DD outlier windows')
                else:
                    df_dxydd.to_csv(file_basename+'DxyDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dxydd.to_csv(file_basename+'DxyDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # select Fst and DD double outliers 
                df_fstdd = df_dat[(df_dat.fstdd == 2)]
                if df_fstdd.empty:
                    print('No Fst DD double outlier windows')
                else:
                    df_fstdd.to_csv(file_basename+'FstDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_fstdd.to_csv(file_basename+'FstDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # # select AFD_1_2, Fst, DD triple outliers 
                # df_afddxyfst = df_dat[(df_dat.afddxyfst == 3)]
                # if df_afddxyfst.empty:
                #     print('No AFD('+args.coh1+'-'+args.coh2+') Dxy Fst triple outlier windows')
                # else:
                #     df_afddxyfst.to_csv(file_basename+'afd12DxyFst_outliers'+args.suf+'.csv', index=False)
                #     # write bedfile for gene retrieval and orthologous gene onthology match
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afddxyfst.to_csv(file_basename+'afd12DxyFst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # # select AFD_1_2, Dxy, DD triple outliers 
                # df_afddxydd = df_dat[(df_dat.afddxydd == 3)]
                # if df_afddxydd.empty:
                #     print('No AFD('+args.coh1+'-'+args.coh2+') Dxy DD triple outlier windows')
                # else:
                #     df_afddxydd.to_csv(file_basename+'afd12DxyDD_outliers'+args.suf+'.csv', index=False)
                #     # write bedfile for gene retrieval and orthologous gene onthology match
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afddxydd.to_csv(file_basename+'afd12DxyDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                    
                # # select AFD_1_2, Fst, DD triple outliers 
                # df_afdfstdd = df_dat[(df_dat.afdfstdd == 3)]
                # if df_afdfstdd.empty:
                #     print('No AFD('+args.coh1+'-'+args.coh2+') Fst DD triple outlier windows')
                # else:
                #     df_afdfstdd.to_csv(file_basename+'afd12FstDD_outliers'+args.suf+'.csv', index=False)
                #     # write bedfile for gene retrieval and orthologous gene onthology match
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afdfstdd.to_csv(file_basename+'afd12FstDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # # select AFD_2_1, Dxy, DD triple outliers 
                # df_afd21dxyfst = df_dat[(df_dat.afd21dxyfst == 3)]
                # if df_afd21dxyfst.empty:
                #     print('No AFD('+args.coh2+'-'+args.coh1+') Dxy Fst triple outlier windows')
                # else:
                #     df_afd21dxyfst.to_csv(file_basename+'afd21DxyFst_outliers'+args.suf+'.csv', index=False)
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afd21dxyfst.to_csv(file_basename+'afd21DxyFst_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
                                        
                # select Dxy, Fst, DD triple outliers 
                df_dxyfstdd = df_dat[(df_dat.dxyfstdd == 3)]
                if df_dxyfstdd.empty:
                    print('No Dxy Fst DD triple outlier windows')
                else:
                    df_dxyfstdd.to_csv(file_basename+'DxyFstDD_outliers'+args.suf+'.csv', index=False)
                    # write bedfile for gene retrieval and orthologous gene onthology match
                    header = ["scaffold", "bedstart", "end"]
                    df_dxyfstdd.to_csv(file_basename+'DxyFstDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)

                # # select AFD12, Dxy, Fst, DD quadruple outliers 
                # df_afddxyfstdd = df_dat[(df_dat.afddxyfstdd == 4)]
                # if df_afddxyfstdd.empty:
                #     print('No AFD('+args.coh1+'-'+args.coh2+') Dxy Fst DD quadruple outlier windows')
                # else:
                #     df_afddxyfstdd.to_csv(file_basename+'afd12DxyFstDD_outliers'+args.suf+'.csv', index=False)
                #     # write bedfile for gene retrieval and orthologous gene onthology match
                #     header = ["scaffold", "bedstart", "end"]
                #     df_afddxyfstdd.to_csv(file_basename+'afd12DxyFstDD_outliers'+args.suf+'.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)


        # set outliers to be inherited to step 6 to generate graphs
        self.out_afd_up = out_afd_up
        self.out_afd_low = out_afd_low
        self.out_dxy = out_dxy
        self.out_fst = out_fst
        self.out_dd = out_dd
        self.c1_c2 = c1_c2


    def step6(self):
        ###### STEP 6 ######
        ## create histograms with the CI cutoff value for Dxy, Fst, DD and window length
        args = self.args
        if args.o is 'na':
            outputdir = str(args.coh1+args.coh2)
        else:
            outputdir = str(args.o+args.coh1+args.coh2)

        contrast = args.coh1 + args.coh2
        args = self.args
        self.outputdir = outputdir
        self.contrast = contrast
        cwd = self.cwd
        print('\n\tSTEP 6: Create population genetics histograms')
        graphs = []
        for dirName, subdirList, fileList in os.walk(outputdir):
            for file in fileList:
                if file.startswith(contrast) and file.endswith(str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_allsites'+args.suf+'.csv') == True:
                    graphs.append(dirName+'/'+file)
        graphs_sorted = natsorted(graphs)

        if os.path.exists(outputdir+'/graphs') == False:
            os.mkdir(outputdir+'/graphs')

        print('Creating histograms for Dxy, Fst, DD and window length')
        for file in graphs_sorted:
            count = 0
            rfile = open(outputdir+'/graphs/'+contrast+'_graphs.r', 'w')
            rfile.write('# load table\n'+
                        'test <- read.csv("'+cwd+'/'+file+'", header=T)\n'+
                        '#str(test)\n'+
                        '#load library\n'+
                        'library(ggplot2)\n'+
                        'library(methods)\n'+
                        'p.afd <- qplot('+str(self.c1_c2)+', data=test, binwidth=0.001) + geom_vline(xintercept='+str(self.out_afd_up)+', color="grey30", linetype="dashed") + geom_vline(xintercept='+str(self.out_afd_low)+', color="grey30", linetype="dotted") + theme_bw()\n'+
                        'p.dxy <- qplot(Dxy, data=test, binwidth=0.001, xlim=c(0, (max(test$Dxy)+0.01))) + geom_vline(xintercept='+str(self.out_dxy)+', color="grey30", linetype="dashed") + theme_bw()\n'+
                        'p.fst <- qplot(Fst, data=test, binwidth=0.001, xlim=c(0, (max(test$Fst)+0.01))) + geom_vline(xintercept='+str(self.out_fst)+', color="grey30", linetype=2) + theme_bw()\n'+
                        'p.dd <- qplot(DD, data=test, binwidth=0.001) + geom_vline(xintercept='+str(self.out_dd)+', color="grey30", linetype=2) + theme_bw()\n'+
                        'p.length <- qplot(length_bp, data=test, binwidth=100, xlim=c(0,26560)) + theme_bw()\n\n'+
                        'p.pi_'+args.coh1+' <- qplot('+args.coh1+'_pi, data=test, binwidth=0.001) + theme_bw()\n\n'+
                        'p.pi_'+args.coh2+' <- qplot('+args.coh2+'_pi, data=test, binwidth=0.001) + theme_bw()\n\n'+
                        'p.pi_afd <- qplot(absdiff, '+args.coh1+'_pi, data=test, xlim=c(0,1)) + theme_bw() + geom_smooth(method="lm")\n\n'+
                        '# export as pdf\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_AFD_histogram'+args.suf+'.pdf")\n'+
                        'p.afd\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_Dxy_histogram'+args.suf+'.pdf")\n'+
                        'p.dxy\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_Fst_histogram'+args.suf+'.pdf")\n'+
                        'p.fst\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_DD_histogram'+args.suf+'.pdf")\n'+
                        'p.dd\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_length_histogram.pdf")\n'+
                        'p.length\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+args.coh1+'_pi_histogram.pdf")\n'+
                        'p.pi_'+args.coh1+'\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+args.coh2+'_pi_histogram.pdf")\n'+
                        'p.pi_'+args.coh2+'\n'+
                        'dev.off()\n'+
                        'pdf(file="'+outputdir+'/graphs/'+contrast+'_'+str(args.snps)+'SNPs_'+args.coh1+'_pi_afd.pdf")\n'+
                        'p.pi_afd\n'+
                        'dev.off()')
            rfile.close()

            cmd = ('Rscript '+outputdir+'/graphs/'+args.coh1+args.coh2+'_graphs.r')
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

            count += 1
            os.remove(outputdir+'/graphs/'+args.coh1+args.coh2+'_graphs.r')
        # print('\n\tDONE\n')    
    

    # define the sequence of steps that should be run
    def step1to6(self):
        self.step1()
        self.step2()
        self.step3()
        self.step4()
        self.step5()
        self.step6()
    def step4to6(self):
        self.step4()
        self.step5()
        self.step6()
    def step5to6(self):
        self.step5()
        self.step6()

if __name__ == '__main__':
    import os, sys