## genome scan pipeline, written by Jeff DaCosta, streamlined by Christian Sailer
## September 2016, updated 13 February 2017

import os, sys, argparse
from PerPopMetrics_classes_1 import G1
# from PerPopMetrics_classes_interval import Ginterval

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Estimation of per population population gentic metrics nucleotide diversity (π), Wattersons theta, Tajimas D '+
                                             'and Fai and Wus H. The input files are joint-genotyped and best practise filtered allele count tables.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full or relative path to input AC table directory')
parser.add_argument('-interval', type=str, metavar='interval_bed', default='na', help='Full or relative path to interval bed file')
parser.add_argument('-snps', type=int, metavar='snps_per_window', default='100', help='Number of variant sites (SNPs) per window [100]')
parser.add_argument('-per', type=float, metavar='percent_missing_data', default='0.2', help='Ratio of missing data allowed, eg 0.2 allows 80percent missing data [0.25]')
parser.add_argument('-o', type=str, metavar='outputdir_path', default='na', help='Optional output directory relative path, contrast directory will be a subdir of this one')
parser.add_argument('-suf', type=str, metavar='suffix_species', default='', help='Suffix to append to the file name, 2-letter species abbreviation, eg. _Aa')
args = parser.parse_args()

print('\n########################################################\n\t\tDivergenceScan pipeline\n\n\t\tDESCRIBING POPULATIONS\n\n#########################################################')

if args.interval is 'na':
    if args.o is 'na':
        outputdir = str('PerPopMetrics_Result/')
    else:
        outputdir = str(args.o)
else:
    # obtain interval naming
    tspace = args.interval.split(sep='/')
    space = tspace[len(tspace)-1]
    # define outputdir
    if args.o is 'na':
        outputdir = str('PerPopMetrics_Result_'+space+'/')
    else:
        outputdir = str(args.o+'_'+space+'/')


# if args.o is 'na':
#     outputdir = str('PerPopMetrics_Result/')
# else:
#     outputdir = str(args.o+'PerPopMetrics_Result/')
#     if os.path.exists(args.o) == False:
#         os.mkdir(args.o)

# # check if output directory exists, create it if necessary
# if os.path.exists(outputdir) == False:
#     os.mkdir(outputdir)
#     print('\n\tCreated directory '+outputdir)

# which input populations?
ACs = []
for dirName, subdirList, fileList in os.walk(args.i):
    for file in fileList:
        if file.endswith('_raw.table') == True:
            ACs.append(dirName+'/'+file)

pops = []
for i in range(len(ACs)):
    tin = ACs[i].split(sep='/')
    inname = tin[len(tin)-1].split(sep='_')[0]
    pops.append(inname)
lastpop = pops[len(ACs)-1]

###### EXECUTION ######
if __name__ == '__main__':
    

    # Two main analysis streams based on set interval file
    if args.interval is 'na':
        genome_scan = G1(args=args)
        if os.path.isfile(outputdir+'/'+lastpop+'_AC'+args.suf+'.table') == False:
            genome_scan.step1to4()
        elif os.path.isfile(outputdir+'/'+lastpop+'_metrics_'+str(args.snps)+'SNPs'+args.suf+'.txt') == False:
            genome_scan.step3to4()
        elif os.path.exists(outputdir+'/graphs/'+inname+'_length_bp.pdf') == False:
            genome_scan.step4to4()


    else:
        genome_scan = Ginterval(args=args)
        if os.path.isfile(outputdir+'/'+lastpop+'_AC'+args.suf+'.table') == False:
            genome_scan.stepsubsetto4()
        elif os.path.isfile(outputdir+'/'+lastpop+'_metrics_'+str(args.snps)+'SNPs'+args.suf+'.txt') == False:
            genome_scan.stepsubset3to4()
        elif os.path.exists(outputdir+'/graphs/'+inname+'_length_bp.pdf') == False:
            genome_scan.step4to4()




       
    # print('\n Succesfully found selective sweep candidates for '+args.coh1+'_'+args.coh2+' contrast and '+str(args.snps)+' SNP windows\n\n')
    # print('\n\tDONE\n')


