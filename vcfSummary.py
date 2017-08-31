# Author: Liping Hou
# Email:  houliping@gmail.com
# March 9th, 2015

import sys
import argparse
import datetime 
import gzip
import math

parser = argparse.ArgumentParser(description="Get variant-level (AC, AN, missing rate, ...) and individual-level (Ti/Tv, number of heterozygotes, ...) summary from a VCF file")
parser.add_argument('input', help="The VCF input file")
parser.add_argument('output', help='The prefix of output files')
#parser.add_argument('--by', help='Get summary results by variant/individual/both (default: both)', choices=['variant','individual','both'], default='both')
args = parser.parse_args()

def vcfSummary(vcf, file):
#def vcfSummary(vcf, file, by):
    print()
    print("Analysis started: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()
    print("Reading vcf from [ {} ]".format(vcf))
    if vcf.endswith(".gz"):
        input=gzip.open(vcf,'rt')
    else:
        input = open(vcf)
    output_indv = open(file+'.ind', 'w')
    output_var = open(file+'.var', 'w')
    iMiss = {}
    Ti = {}
    Tv = {}
    NALT = {}
    NHET = {}
    NVAR = {}
    DEPTH = {}
    QUAL = {}
    AD_REF = {}
    AD_ALT = {}
#To do: add codes to estimate allele balance of all HET at the individual level
    def TiTv(x, y):
        nucleotides = ['A', 'C', 'T', 'G']
        if x not in nucleotides:
            print("{} is not a 'A/T/C/G'".format(x))
            sys.exit()
        if y not in nucleotides:
            print("{} is not a 'A/T/C/G'".format(y))
            sys.exit()
        if (x == 'A' and y == 'G') or (x == 'G' and y == 'A') or (x == 'T' and y == 'C') or (x == 'C' and y == 'T'):
            return 'Ti'
        else:
            return 'Tv'

    def genoParser(data, format):
        f = format.split(':')
        g = data.split(':')
        geno = {}
        for i in range(len(f)):
            geno[f[i]] = g[i]
        return(geno)

    nMarker = nMAP = nMNP = nIns = nDel = nSNP = nTi = nTv = 0
    
    output_var.write("CHROM\tPOS\tID\tREF\tALT\tAC\tAF\tAN\tMISS\tAA/AB/BB/NN\n")
    
    for line in input:
        if line.startswith("##"): # skip vcf Meta-information lines
            continue
        elif line.startswith("#"): 
            header = line[1:].strip()
            ids = header.split()[9:]
            nIndiv = len(ids)
            print("Found {} individuals from [ {} ]".format(nIndiv, vcf))
            print("Writing variant level summary results to: [ {} ]".format(file + '.var'))
            for i in ids:
                iMiss[i] = 0
                Ti[i] = 0
                Tv[i] = 0
                NALT[i] = 0
                NHET[i] = 0
                NVAR[i] = 0
                DEPTH[i] = 0
                QUAL[i] = 0
                AD_REF[i] = 0
                AD_ALT[i] = 0
        else: 
            data = line.strip().split()
            geno_format = data[8]
            nMarker += 1
            if len(data[4].split(",")) != 1:
                nMAP += 1
            else:
                if len(data[3]) > len(data[4]):
                    nDel += 1
                elif len(data[3]) < len(data[4]):
                    nIns += 1
                elif len(data[3]) == len(data[4]) and len(data[3]) > 1:
                    nMNP += 1
                elif len(data[3]) == 1 and len(data[4]) == 1:
                    nSNP += 1
                    TiTv_status = TiTv(data[3], data[4])
                    if TiTv_status == 'Ti':
                        nTi += 1
                    else:
                        nTv += 1
                A0 = A1 = A2 = 0
                AA = AB = BB = NN = 0
                AC = AF = AN = MISS = 0 # Have to initialize these variables or will cause trouble when all genotypes of one variant start with ./.
                if len(data[9].split(":")[0]) == 3: #if diploid then continue else skip
                    if len(data[3]) == 1 and len(data[4]) == 1: # If 'SNP' then estimate Ti/Tv and iMiss, NVAR, NALT, NHET
                        TiTv_status = TiTv(data[3], data[4])
                        for i in range(9,len(data)):
                            if data[i].startswith("./."):
                                A0 += 2
                                NN += 1
                                iMiss[ids[i-9]] += 1
                            elif data[i][0] == '.' or data[i][2] == '.':
                                A0 += 2
                                NN += 1 # count half-calls as NNs
                                iMiss[ids[i-9]] += 1
                            elif data[i][0] == '0' and data[i][2] == '0':
                                A1 += 2
                                AA += 1
                                NVAR[ids[i-9]] += 1
                                #data[i] = "0/0:.:.:.:.:.,."
                            else:
                                #print(data[i])
                                geno = genoParser(data[i], geno_format)
                                #GT:FT:GQ:GL:DP:AD
                                GT = geno['GT']
                                try:
                                    DP = int(geno['DP'])
                                except ValueError:
                                    DP = 0
                                try:
                                    GQ = float(geno['GQ'])
                                except ValueError:
                                    GQ = 0.0
                                if GT.split(GT[1])[0] != GT.split(GT[1])[1]:
                                    A1 += 1
                                    A2 += 1
                                    AB += 1
                                    NVAR[ids[i-9]] += 1
                                    NHET[ids[i-9]] += 1
                                    NALT[ids[i-9]] += 1
                                    DEPTH[ids[i-9]] += DP
                                    QUAL[ids[i-9]] += GQ
                                    if TiTv_status == 'Ti':
                                        Ti[ids[i-9]] += 1
                                    else:
                                        Tv[ids[i-9]] += 1
                                #elif GT == '1/1' or GT == '1|1':
                                elif GT.split(GT[1])[0] == GT.split(GT[1])[1] and GT.split(GT[1])[0] == '1':
                                    A2 += 2
                                    BB += 1
                                    NVAR[ids[i-9]] += 1
                                    NALT[ids[i-9]] += 1
                                    DEPTH[ids[i-9]] += DP
                                    QUAL[ids[i-9]] += GQ
                                    if TiTv_status == 'Ti':
                                        Ti[ids[i-9]] += 1
                                    else:
                                        Tv[ids[i-9]] += 1
                    else: # if 'non-SNP' then only estimate iMiss, NVAR, NALT, NHET
                        for i in range(9,len(data)):
                            if data[i].startswith("./."):
                                A0 += 2
                                NN += 1
                                iMiss[ids[i-9]] += 1
                            elif data[i][0] == '.' or data[i][2] == '.':
                                A0 += 2
                                NN += 1 # count half-calls as NNs
                                iMiss[ids[i-9]] += 1
                            elif data[i][0] == '0' and data[i][2] == '0':
                                A1 += 2
                                AA += 1
                                NVAR[ids[i-9]] += 1
                                #data[i] = "0/0:.:.:.:.:.,."
                            else:
                                #print(data[i])
                                geno = genoParser(data[i], geno_format)
                                #GT:FT:GQ:GL:DP:AD
                                GT = geno['GT']
                                try:
                                    DP = int(geno['DP'])
                                except ValueError:
                                    DP = 0
                                try:
                                    GQ = float(geno['GQ'])
                                except ValueError:
                                    GQ = 0.0
                                #DP = int(geno['DP'])
                                #GQ = float(geno['GQ'])
                                if GT.split(GT[1])[0] != GT.split(GT[1])[1]:
                                    A1 += 1
                                    A2 += 1
                                    AB += 1
                                    NVAR[ids[i-9]] += 1
                                    NHET[ids[i-9]] += 1
                                    NALT[ids[i-9]] += 1
                                    DEPTH[ids[i-9]] += DP
                                    QUAL[ids[i-9]] += GQ
                                #elif GT == '1/1' or GT == '1|1':
                                elif GT.split(GT[1])[0] == GT.split(GT[1])[1] and GT.split(GT[1])[0] == '1':
                                    A2 += 2
                                    BB += 1
                                    NVAR[ids[i-9]] += 1
                                    NALT[ids[i-9]] += 1
                                    DEPTH[ids[i-9]] += DP
                                    QUAL[ids[i-9]] += GQ
                AC = A2
                AN = A1 + A2
                if AN != 0:
                    AF = A2/AN
                else:
                    AF = math.inf 
                MISS = NN / nIndiv
                output_var.write("{}\t{}\t{:.3f}\t{}\t{:.3f}\t{}/{}/{}/{}\n".format("\t".join(str(j) for j in data[0:5]), AC, AF, AN, MISS, AA, AB, BB, NN))
    print("Writing individual level summary results to: [ {} ]".format(file + '.ind'))
    output_indv.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('ID', 'NVAR', 'NALT', 'NHET', 'MISS', 'TiTv','QUAL','DP'))
    for i in ids:
        output_indv.write("{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.2f}\t{:.2f}\n".format(i, NVAR[i], NALT[i], NHET[i], iMiss[i]/nMarker, Ti[i]/Tv[i], QUAL[i]/NALT[i], DEPTH[i]/NALT[i]))
    input.close()
    output_var.close()
    output_indv.close()
    print("Found {} markers from [ {} ]".format(nMarker, vcf))
    print("Including {} multi-allelic markers, {} MNPs, {} Insertions, {} Deletions, {} SNPs".format(nMAP, nMNP, nIns, nDel, nSNP))
    print("The overall Ti/Tv ratio is {:.3f}".format(nTi/nTv))
    print()
    print("Analysis finished: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()

vcfSummary(args.input, args.output)
#vcfSummary(args.input, args.output, args.by)
