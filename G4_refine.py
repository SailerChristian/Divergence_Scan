## genome scan pipeline, refining
## by Christian Sailer, 23 July 2017

import os, sys, argparse
from GS_classes_2 import G2

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Refining the candidate gene lists from steps GS1 to GS2. '+
                                 'Output are the metrics per candiate gene plus SNPeff annotations per gene.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full or relative path to input AC table directory')
parser.add_argument('-coh1', type=str, metavar='cohort_1', required=True, help='REQUIRED: Name of first cohort')
parser.add_argument('-coh2', type=str, metavar='cohort_2', required=True, help='REQUIRED: Name of second cohort')
parser.add_argument('-snps', type=int, metavar='snps_per_window', default='100', help='Number of SNPs per window [100]')
parser.add_argument('-o', type=str, metavar='outputdir_path', default='na', help='Optional output directory relative path, contrast directory will be a subdir of this one')
parser.add_argument('-cut', type=float, metavar='cutoff_ratio', default='0.5', help='Top percentile, defines outliers [0.5]')
parser.add_argument('-suf', type=str, metavar='suffix_species', default='', help='Suffix to append to the file name, 2-letter species abbreviation, eg. _Aa')
parser.add_argument('-ovlp', type=float, metavar='percent_overlap', default='0.0001', help='Precentage of base pairs that overlap the search pattern as percentage [0.0001]')
parser.add_argument('-outlier', type=str, metavar='outliers_to_plot', default='na', help='Define for which outlier combinations to produce graphs, eg DxyDD [na]') 
args = parser.parse_args()

#cwd = os.getcwd()
#Test whether specified -o exist, create if necessary
if args.o is 'na':
    outputdir = str(args.coh1+args.coh2+'/refined/')
else:
    outputdir = str(args.o+args.coh1+args.coh2+'/refined/')
    if os.path.exists(outputdir) == False:
        os.mkdir(outputdir)

# Test whether sub directory of outputdir exists, create if it does not
if os.path.exists(outputdir+args.outlier+'/') == False:
    os.mkdir(outputdir+args.outlier+'/')
    print('\nCreated directory '+outputdir+args.outlier+'/')
contrast = args.coh1+args.coh2

###### EXECUTION ######
if __name__ == '__main__':
    GS = G2(args=args)
    if os.path.isfile(outputdir+args.outlier+'/'+contrast+'_'+str(args.snps)+'SNPs_'+str(int(10000*args.cut))+'ppm_'+str(args.outlier)+'_'+str(int(args.ovlp))+'ol_refined_AC.table') == False:
        GS.step1to5()
    # else:
    #     print('Jump')
    #     GS.step3_gene()

