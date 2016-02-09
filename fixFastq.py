# Author: Liping Hou
# Email:  liping.hou@nih.gov
# Feb, 2016

import sys
import gzip 
import argparse
import datetime 

parser = argparse.ArgumentParser(description="Fix FASTQ files (remove singletons, resolve pairs, recode quality scores to Illumina-1.8 if needed)")
parser.add_argument('fastq1', help="input fastq file for first read in paired end data")
parser.add_argument('fastq2', help='input fastq file for second read in paired end data')
parser.add_argument('output', help='the prefix of output fastq files')
parser.add_argument('--checkEncodingOnly', help='Use this flag if you only want to check the encoding', action='store_true')
args = parser.parse_args()

def fixFastq(fastq1, fastq2, out, encodingOnly):
    print()
    print("Analysis started: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()
    if fastq1.endswith(".gz"):
        f1 = gzip.open(fastq1, 'rt')
    else:
        f1 = open(fastq1)
    if fastq2.endswith(".gz"):
        f2 = gzip.open(fastq2, 'rt')
    else:
        f2 = open(fastq2)
    try:
        encodingOnly
    except NameError:
        encodingOnly = False
    ReadNames1 = []
    SEQ1 = {}
    QUAL1 = {}
    ReadNames2 = []
    SEQ2 = {}
    QUAL2 = {}
    MinQual1 = MinQual2 = 126
    MaxQual1 = MaxQual2 = 33

    #Notes: https://en.wikipedia.org/wiki/FASTQ_format
    #    'Illumina-1.3': [64, 104]
    #    'Illumina-1.5': [67, 104]
    #    'Illumina-1.8': [33, 93]
    
    def getEncoding(qual_min, qual_max):
        if qual_min >= 33 and qual_max <= 93:
            return "Illumina-1.8"
        elif qual_min >= 66 and qual_max <= 105:
            return "Illumina-1.5"
        elif qual_min >= 64 and qual_max <= 104:
            return "Illumina-1.3"
        else:
            return "Not a valid Illumina encoding!"

    print("Reading the first fastq file from [ {} ]".format(fastq1))
    if encodingOnly:
        for line in f1:
            res = line.strip()
            if res.startswith("@"):
                ReadNames1.append(res)
                next(f1)
                next(f1)
                qual1 = next(f1).strip()
                vals = [ord(m) for m in qual1]
                MinQual1 = min(vals) if min(vals) < MinQual1 else MinQual1
                MaxQual1 = max(vals) if max(vals) > MaxQual1 else MaxQual1
    else:
        for line in f1:
            res = line.strip()
            if res.startswith("@"):
                ReadNames1.append(res)
                SEQ1[res] = next(f1)
                next(f1)
                qual1 = next(f1).strip()
                QUAL1[res] = qual1
                vals = [ord(m) for m in qual1]
                MinQual1 = min(vals) if min(vals) < MinQual1 else MinQual1
                MaxQual1 = max(vals) if max(vals) > MaxQual1 else MaxQual1
    R1 = set(ReadNames1)
    print("Processed {} reads (including {} unique reads)".format(len(ReadNames1), len(R1)))
    print("The range of quality scores is [{}, {}] (Guessed encoding: {})".format(MinQual1,MaxQual1,getEncoding(MinQual1, MaxQual1)))
    print("Reading the second fastq file from [ {} ]".format(fastq2))
    if encodingOnly:
        for line in f2:
            res = line.strip()
            if res.startswith("@"):
                ReadNames2.append(res)
                next(f2)
                next(f2)
                qual2 = next(f2).strip()
                vals = [ord(m) for m in qual2]
                MinQual2 = min(vals) if min(vals) < MinQual2 else MinQual2
                MaxQual2 = max(vals) if max(vals) > MaxQual2 else MaxQual2
    else:
        for line in f2:
            res = line.strip()
            if res.startswith("@"):
                ReadNames2.append(res)
                SEQ2[res] = next(f2)
                next(f2)
                qual2 = next(f2).strip()
                QUAL2[res] = qual2
                vals = [ord(m) for m in qual2]
                MinQual2 = min(vals) if min(vals) < MinQual2 else MinQual2
                MaxQual2 = max(vals) if max(vals) > MaxQual2 else MaxQual2
    R2 = set(ReadNames2)
    print("Processed {} reads (including {} unique reads)".format(len(ReadNames2), len(R2)))
    print("The range of quality scores is [{}, {}] (Guessed encoding: {})".format(MinQual2,MaxQual2,getEncoding(MinQual2, MaxQual2)))
    if not encodingOnly:
        paired = set.intersection(R1, R2)
        print("Found {} read pairs from {} and {}".format(len(paired), fastq1, fastq2))

    def q33(qual):
        vals = [chr(ord(m)-31) for m in qual]
        return ''.join(n for n in vals)

    if encodingOnly:
        f1.seek(0,0)
        f2.seek(0,0)
        if getEncoding(MinQual1,MaxQual1) == 'Illumina-1.8': # Assume the encoding is the same for FASTQ2
            print("No need to fix the fastq files")
        elif getEncoding(MinQual1,MaxQual1) == 'Illumina-1.3' or getEncoding(MinQual1,MaxQual1) == 'Illumina-1.5':
            output1 = open(out+'_R1.fastq', 'w')
            output2 = open(out+'_R2.fastq', 'w')
            print("Writing fixed fastq files to {} and {}".format(out+'_R1.fastq', out+'_R2.fastq'))
            for line in f1:
                rid = line
                seq=next(f1)
                line3=next(f1)
                qual=next(f1).strip()
                output1.write("{}".format(rid))
                output1.write("{}".format(seq))
                output1.write("{}".format(line3))
                output1.write("{}\n".format(q33(qual)))
            for line in f2:
                rid = line
                seq=next(f2)
                line3=next(f2)
                qual=next(f2).strip()
                output2.write("{}".format(rid))
                output2.write("{}".format(seq))
                output2.write("{}".format(line3))
                output2.write("{}\n".format(q33(qual)))
            output1.close()
            output2.close()
        else:
            print("Error: invalid encoding of the FASTQ files!")
    else:
        output1 = open(out+'_R1.fastq', 'w')
        output2 = open(out+'_R2.fastq', 'w')
        print("Writing fixed fastq files to {} and {}".format(out+'_R1.fastq', out+'_R2.fastq'))
        if getEncoding(MinQual1,MaxQual1) == 'Illumina-1.8':
            for i in paired:
                output1.write("{}\n".format(i))
                output1.write("{}".format(SEQ1[i]))
                output1.write("+\n")
                output1.write("{}\n".format(QUAL1[i]))
                output2.write("{}\n".format(i))
                output2.write("{}".format(SEQ2[i]))
                output2.write("+\n")
                output2.write("{}\n".format(QUAL2[i]))
        elif getEncoding(MinQual1,MaxQual1) == 'Illumina-1.3' or getEncoding(MinQual1,MaxQual1) == 'Illumina-1.5':
            for i in paired:
                output1.write("{}\n".format(i))
                output1.write("{}".format(SEQ1[i]))
                output1.write("+\n")
                output1.write("{}\n".format(q33(QUAL1[i])))
                output2.write("{}\n".format(i))
                output2.write("{}".format(SEQ2[i]))
                output2.write("+\n")
                output2.write("{}\n".format(q33(QUAL2[i])))
        else:
            print("Error: invalid encoding of the FASTQ files!")
        output1.close()
        output2.close()
    f1.close()
    f2.close()
    print()
    print("Analysis finished: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()

fixFastq(args.fastq1, args.fastq2, args.output, args.checkEncodingOnly)
