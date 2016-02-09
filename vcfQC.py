# Author: Liping Hou
# Email:  houliping@gmail.com
# March 4th, 2015

import argparse
import datetime 
import sys
import gzip

parser = argparse.ArgumentParser(description="QC the vcf file by: 1) setting genotypes with GQ/DP/AB below the specified thresholds to missing; 2) removing/marking variants with high missing rate or low allele count")
parser.add_argument('input', help="Name of the input file (VCF format)")
parser.add_argument('output', help='Name of the output file with all bad genotypes excluded')
parser.add_argument('--caller', help='which variant caller was used to generate the VCF input file', choices=['gatk','freebayes','platypus','samtools'], default='gatk')
parser.add_argument('--GQ', help='Genotypes with "Genotype Quality" < N will be excluded (default: 20)', metavar="N", type=int, default=20)
parser.add_argument('--DP', help='Genotypes with "Total Depth" < N will be excluded (default: 8)', metavar="N", type=int, default=8)
parser.add_argument('--AB', help='Heterozygotes with "Allele Balance Ratio" < N (or > 1 - N) will be excluded (default: 0.25)', metavar="N", type=float, default=0.25)
parser.add_argument('--MISSING', help='Add a "HighMissing" flag to the FILTER column for all Variants with missing rate > N (default: 0.25)', metavar="N", type=float, default=0.25)
parser.add_argument('--AC', help='Add a "LowAC" flag to the FILTER column for all variants with AC < N (default: 1)', metavar="N", type=int, default=1)
parser.add_argument('--RemoveFiltered', help='Use this flag to remove all "HighMissing/LowAC" variants', action='store_true')
args = parser.parse_args()

def vcfQC(vcf, file, caller, GQ_threshold, DP_threshold, AB_threshold, missing_rate, allele_count, removeFiltered):
    print()
    print("Analysis started: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()
    print("Reading vcf from [ {} ]".format(vcf))
    if vcf.endswith(".vcf.gz"):
        input = gzip.open(vcf, 'rt')
    else:
        input = open(vcf)
    output = open(file, 'w')

    try:
        removeFiltered
    except NameError:
        removeFiltered = False

    def genoParser(data, format):
        f = format.split(':')
        g = data.split(':')
        geno = {}
        for i in range(len(f)):
            geno[f[i]] = g[i]
        return(geno)

    badGQ = badDP = badAB = 0
    nMarkers = nHalfCalls = badGeno = highMiss = lowAC = nMAP = nHet = nHaploid = nHetbadAD = nHetbadGQ = nHetbadDP = 0
    for line in input:
        if line.startswith("##"): # write meta-information lines to output 
            output.write(line)
        elif line.startswith("#CHROM"):
            if not removeFiltered:
                output.write("##FILTER=<ID=HighMissing,Description=\"Missing rate is greater than {}\">\n".format(missing_rate))
                output.write("##FILTER=<ID=LowAC,Description=\"AC is less than {}\">\n".format(allele_count))
            nIndiv = len(line.strip().split()) - 9
            output.write(line)
        else:
            data = line.strip().split()
            geno_format = data[8]
            nMarkers += 1
            if len(data[4].split(",")) != 1: # remove multi-allelic markers (MAPs)
                nMAP += 1
                continue
            A2 = 0 #alternative allele counter
            NN = 0 #missing genotype counter for each variant
            if len(data[9].split(":")[0]) == 3: # diploid only!!!
                for i in range(9,len(data)):
                    if data[i][0] == '.' and data[i][2] == '.': # if noCall
                        badGeno += 1
                        NN += 1
                    elif data[i][0] == '.' or data[i][2] == '.': # if halfCall
                        badGeno += 1
                        NN += 1 # count half-calls as NNs
                        nHalfCalls += 1
                        data[i] = './.' + data[i][3:] # set all half-calls to missing
                    elif data[i][0] == data[i][2]: # if hom
                        #print(data[i])
                        geno = genoParser(data[i], geno_format)
                        GT = geno['GT']
                        if geno['GQ'] == '.':
                            GQ = None
                        else:
                            GQ = float(geno['GQ'])
                        if caller == 'gatk':
                            if geno['DP'] == ".":
                                DP = None
                            else:
                                DP = int(geno['DP'])
                        elif caller == 'freebayes':
                            if geno['DP'] == ".":
                                DP = None
                            else:
                                DP = int(geno['DP'])
                        elif caller == 'platypus':
                            if geno['NR'] == ".":
                                DP = None
                            else:
                                DP = int(geno['NR'])
                        else: # Unfortunately samtools does not provide DP and AD
                            DP = DP_threshold + 1
                        if GQ == None or GQ < GQ_threshold:
                            geno['GT'] = './.'
                            badGQ += 1
                            NN += 1
                        elif DP == None or DP < DP_threshold:
                            geno['GT'] = './.'
                            badDP += 1
                            NN += 1
                        data[i] = ":".join(geno[i] for i in geno_format.split(":"))
                        GT = geno['GT']
                        #if GT == '1/1' or GT == '1|1':
                        if GT.split(GT[1])[0] == GT.split(GT[1])[1] and GT.split(GT[1])[0] == '1':
                            A2 += 2
                    else: #if het
                        #print(data[i])
                        nHet += 1
                        geno = genoParser(data[i], geno_format)
                        GT = geno['GT']
                        if geno['GQ'] == '.':
                            GQ = None
                        else:
                            GQ = float(geno['GQ'])
                        if caller == 'gatk':
                            if geno['DP'] == ".":
                                DP = None
                            else:
                                DP = int(geno['DP'])
                            if len(geno['AD'].split(',')) == 2:
                                AD1 = geno['AD'].split(',')[0]
                                AD2 = geno['AD'].split(',')[1]
                                if AD1 == '.' or AD2 == '.':
                                    AB = None
                                elif int(AD1) + int(AD2) == 0:
                                    AB = None
                                else:
                                    AB = int(AD2)/(int(AD1) + int(AD2))
                            else: # ignore if length of AD > 2
                                AB = 0.5
                                nHetbadAD += 1
                        elif caller == 'freebayes':
                            if geno['DP'] == ".":
                                DP = None
                            else:
                                DP = int(geno['DP'])
                            AD1 = geno['RO']
                            AD2 = geno['AO']
                            if AD1 == '.' or AD2 == '.':
                                AB = None
                            elif int(AD1) + int(AD2) == 0:
                                AB = None
                            else:
                                AB = int(AD2)/(int(AD1) + int(AD2))
                        elif caller == 'platypus':
                            if geno['NR'] == ".":
                                DP = None
                            else:
                                DP = int(geno['NR'])
                            AD2 = geno['NV']
                            if DP == None or AD2 == '.':
                                AB = None
                            AB = int(AD2)/int(DP)
                        else: # Unfortunately samtools does not provide DP and AD
                            DP = DP_threshold + 1
                            AB = 0.5
                        if GQ == None or GQ < GQ_threshold:
                            geno['GT'] = './.'
                            badGQ += 1
                            nHetbadGQ += 1
                            NN += 1
                        elif DP == None or DP < DP_threshold:
                            geno['GT'] = './.'
                            badDP += 1
                            nHetbadDP += 1
                            NN += 1
                        elif AB == None or AB > (1-AB_threshold) or AB < AB_threshold:
                            #print("{}\t{}\t{}\t{}".format(geno['GT'],geno['AD'],geno['DP'],AB))
                            geno['GT'] = './.'
                            badAB += 1
                            NN += 1
                        #if GT == '0/1' or GT == '0|1' after QC:
                        GT = geno['GT']
                        if GT.split(GT[1])[0] != GT.split(GT[1])[1]:
                            A2 += 1
                        data[i] = ":".join(geno[i] for i in geno_format.split(":"))
                if A2 < allele_count:
                    lowAC += 1
                    if removeFiltered:
                        pass
                    else:
                        if data[6] == '.' or data[6] == 'PASS':
                            data[6] = 'LowAC'
                        elif 'LowAC' in data[6]:
                            pass
                        else:
                            data[6] = data[6] + ';' + 'LowAC'
                        output.write("{}\n".format("\t".join(str(j) for j in data)))
                elif NN / nIndiv > missing_rate:
                    highMiss += 1
                    if removeFiltered:
                        pass
                    else:
                        if data[6] == '.' or data[6] == 'PASS':
                            data[6] = 'HighMissing'
                        elif 'HighMissing' in data[6]:
                            pass
                        else:
                            data[6] = data[6] + ';' + 'HighMissing'
                        output.write("{}\n".format("\t".join(str(j) for j in data)))
                else:
                    if removeFiltered:
                        if data[6] == '.' or data[6] == 'PASS':
                            output.write("{}\n".format("\t".join(str(j) for j in data)))
                    else:
                        if data[6] == '.':
                            data[6] = "PASS"
                        output.write("{}\n".format("\t".join(str(j) for j in data)))
            else: # non-diploid
                nHaploid += 1
    input.close()
    output.close()
    print("{} markers, {} individuals to be included from [ {} ]".format(nMarkers, nIndiv, vcf))
    print("{} multi-allelic variants were found and excluded".format(nMAP))
    print("{} variants with non-diploid genotypes were excluded".format(nHaploid))
    print("Total missing rate before QC is {}".format(round(badGeno/(nMarkers*nIndiv),3)))
    if nHalfCalls > 0:
        print("{} half-calls were set to missing".format(nHalfCalls))
    print("{}({:.2f}%) genotypes (including: {} heterozygotes) removed with '--GQ {}' option".format(badGQ, badGQ/(nMarkers*nIndiv - badGeno)*100, nHetbadGQ, GQ_threshold))
    print("{}({:.2f}%) genotypes (including: {} heterozygotes) removed with '--DP {}' option".format(badDP, badDP/(nMarkers*nIndiv - badGeno)*100, nHetbadDP, DP_threshold))
    print("{}({:.2f}%) heterozygotes (out of {}) removed with '--AB {}' option".format(badAB, badAB/nHet*100, nHet, AB_threshold))
    if nHetbadAD > 0:
        print("Warning: {} of heterozygotes have malformed AD".format(nHetbadAD))
    print("Total missing rate after QC is {:.2f}".format((badGQ + badDP + badAB + badGeno)/(nMarkers*nIndiv)))
    if removeFiltered:
        print("{}({:.2f}%) variants removed with '--AC {}' option".format(lowAC, lowAC/nMarkers*100, allele_count))
        print("{}({:.2f}%) variants removed with '--MISSING {}' option".format(highMiss, highMiss/nMarkers*100, missing_rate))
    else:
        print("{}({:.2f}%) variants marked as 'lowAC' with '--AC {}' option".format(lowAC, lowAC/nMarkers*100, allele_count))
        print("{}({:.2f}%) variants marked as 'HighMissing' with '--MISSING {}' option".format(highMiss, highMiss/nMarkers*100, missing_rate))
    print("Output file was written to: [ {} ]".format(file))
    print()
    print("Analysis finished: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()

vcfQC(args.input, args.output, args.caller, args.GQ, args.DP, args.AB, args.MISSING, args.AC, args.RemoveFiltered)
