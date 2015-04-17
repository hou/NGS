# Author: Liping Hou
# Email:  houliping@gmail.com
# April 16, 2015

import argparse
import datetime 
import sys
import gzip

parser = argparse.ArgumentParser(description="Backfilling is needed when different VCF files were merged. This script will go through the merged VCF file, and backfill the genotypes to '0/0' when all subejcts provided by the user have missing calls")
parser.add_argument('input', help="Name of the input file (VCF format)")
parser.add_argument('output', help='Name of the output file')
parser.add_argument('--subjects', help='A file that includes groups of subjects for whom the backfilling is needed (This file can have multhiple lines, but subjects that belong to the same group have to be on a single line)',metavar="")
args = parser.parse_args()

def backfillVCF(vcf, file, subjects):
    print()
    print("Analysis started: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()
    if vcf.endswith(".vcf.gz"):
        input = gzip.open(vcf, 'rt')
    else:
        input = open(vcf)
    output = open(file, 'w')
    indv = open(subjects)
    nMarkers = group = 0
    SUBJECT = {}
    SUBJECT_index = {}
    GENO = {}
    BACKFILL= {}
    for line in indv:
        tem = line.strip().split()
        SUBJECT[group] = tem
        BACKFILL[group] = 0
        group += 1

    print("Reading subject list from [ {} ]".format(subjects))
    print("Found {} groups of people from the [ {} ]".format(group, subjects))
    for i in SUBJECT:
        print("Group {} has {} sujects".format(i+1, len(SUBJECT[i])))

    print("Backfilling the vcf file: [ {} ]".format(vcf))
    for line in input:
        if line.startswith("##"): # write meta-information lines to output 
            output.write(line)
        elif line.startswith("#CHROM"): 
            ids = line.strip().split()
            nIndiv = len(ids) - 9
            output.write(line)
            for i in range(group):
                SUBJECT_index[i] = [ids.index(j) for j in SUBJECT[i]]
        else:
            data = line.strip().split()
            nMarkers += 1
            format = data[8].split(":")
            for i in range(group):
                GENO[i] = [data[j][:3] for j in SUBJECT_index[i]]
                if all(k == './.' for k in GENO[i]):
                    BACKFILL[i] += 1
                    for j in SUBJECT_index[i]:
                        data[j] = '0/0' + ":."*len(format)
            output.write("{}\n".format("\t".join(str(j) for j in data)))
    input.close()
    output.close()
    indv.close()
    print("{} markers, {} individuals were found from [ {} ]".format(nMarkers, nIndiv, vcf))
    for i in range(group):
        print("Backfilled {} markers for people of group {}".format(BACKFILL[i], i+1))
    print("Output file was written to: [ {} ]".format(file))
    print()
    print("Analysis finished: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()

backfillVCF(args.input, args.output, args.subjects)
