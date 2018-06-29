## genome scan pipeline, per SNP output
## by Christian Sailer, 7 December 2017
## updated 8 March 2018

import os, sys, argparse
from DS_classes_perSNP import G_NullW

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Refining the candidate gene lists from steps GS1 to GS2. '+
                                 'Output are the metrics per candiate gene plus SNPeff annotations per gene.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full or relative path to directory holding the contrast directory and the input WG AC table directory')
parser.add_argument('-coh1', type=str, metavar='cohort_1', required=True, help='REQUIRED: Name of first cohort')
parser.add_argument('-coh2', type=str, metavar='cohort_2', required=True, help='REQUIRED: Name of second cohort')
parser.add_argument('-o', type=str, metavar='outputdir_path', default='na', help='Optional output directory relative path, contrast directory will be a subdir of this one [na=args.i]')
parser.add_argument('-suf', type=str, metavar='suffix_species', default='', help='Suffix to append to the file name, 2-letter species abbreviation, eg. _Aa')
# parser.add_argument('-genes', type=str, metavar='cand_genes_bed', default='Data/annotations/genes_only.bed', help='Absolute or relative path to bed file containing the candidate gene interval [Data/annotations/genes_only.bed]')
parser.add_argument('-gentxt', type=str, metavar='genenames', default='Data/annotations/all_genes.txt', help='Absolute or relative path to genenames text file [Data/annotations/all_genes.txt]') 
parser.add_argument('-outlier', type=str, metavar='outlier_metric', default='ALL', help='Type of outlier for candidate genes, e.g. ALL, DD, FstDD [ALL]')
args = parser.parse_args()

#cwd = os.getcwd()
#Test whether specified -o exist, create if necessary
if args.o is 'na':
    outputdir = str(args.i+args.coh1+args.coh2+'/perSNP/')
    if os.path.exists(outputdir) == False:
    	os.mkdir(outputdir)
else:
    outputdir = str(args.o+args.coh1+args.coh2+'/perSNP/')
    if os.path.exists(args.o) == False:
        os.mkdir(args.o)
    if os.path.exists(args.o+args.coh1+args.coh2) == False:
        os.mkdir(args.o+args.coh1+args.coh2)
    if os.path.exists(outputdir) == False:
        os.mkdir(outputdir)

print('\n\t###############################')
print('\t## You are using the per SNP Fst output script from the Yant lab.')
print('\t## Please cite our works.\n\t###############################\n')

# # Test whether sub directory of outputdir exists, create if it does not
# if os.path.exists(outputdir+args.outlier+'/') == False:
#     os.mkdir(outputdir+args.outlier+'/')
#     print('\nCreated directory '+outputdir+args.outlier+'/')
contrast = args.coh1+args.coh2
bname = contrast+'_WG'


###### EXECUTION ######
if __name__ == '__main__':
    GS = G_NullW(args=args)
    if os.path.isfile(outputdir+contrast+'_snp_metrics_'+bname+args.suf+'.bed') == False:
        GS.step1to4()
    elif os.path.isfile(outputdir+bname+'_'+args.outlier+'_DS_SNPeff'+args.suf+'.txt') == False:
        GS.step3to4()
    else:
    	print(contrast+' has already been processed for '+args.outlier+'\n')


