# Author: Liping Hou
# Email:  houliping@gmail.com
# Nov. 7th, 2014

import argparse
import datetime 

parser = argparse.ArgumentParser(description="Scan the vcf file for Mendelian Errors")
parser.add_argument('vcf', help="The VCF input file")
parser.add_argument('fam', help='The fam file ')
parser.add_argument('prefix', help='The prefix of the output files')
parser.add_argument('--zeroout','-z', action='store_true', help='Create a new vcf file by zeroing out all Mendelian Errors')
parser.add_argument('--me', metavar="N", type=float, default=0.1, help='Mark all variants with > N Mendelian error rate (based on trios) in the new vcf file')

args = parser.parse_args()

def vcfPedcheck(vcf, fam, prefix, zeroout, me):
    print()
    print("Analysis started: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()
    print("Reading pedigree information from [ {} ]".format(fam))
    ped  = open(fam)
    FID = {}
    PID = {}
    MID = {}
    IID = []
    SEX = {}
    founders = []
    trios = []
    fatherKids = []
    motherKids = []
    nFemale = 0
    nMale = 0
    for indiv in ped:
        tem = indiv.strip().split()
        if tem[4] == '1':
            nMale += 1
        elif tem[4] == '2':
            nFemale += 1
        IID.append(tem[1])
        FID[tem[1]] = tem[0]
        PID[tem[1]] = tem[2]
        MID[tem[1]] = tem[3]
        SEX[tem[1]] = tem[4]
        if tem[2] == "0" and tem[3] == "0":
            founders.append(tem[1])
        elif tem[2] == "0" and tem[3] != "0":
            motherKids.append(tem[1])
        elif tem[2] != "0" and tem[3] == "0":
            fatherKids.append(tem[1])
        else:
            trios.append(tem[1])
    print("{} individuals ({} females, {} males) read from [ {} ]".format(len(IID), nFemale, nMale, fam))
    print("Including {} trios and {} duos".format(len(trios),len(fatherKids)+len(motherKids)))
    print("Reading vcf from [ {} ]".format(vcf))
    try:
        zeroout
    except NameError:
        zeroout = False
    input = open(vcf)
    if zeroout:
        output = open(prefix + '.vcf', 'w')
    mendel = open(prefix + '.mendel', 'w')
    lmendel = open(prefix + '.lmendel', 'w')
    imendel = open(prefix + '.imendel', 'w')
    HetHaploid = open(prefix + '.hh', 'w')
    mendel.write("FID\tPID\tMID\tIID\tCHR\tPOS\tSNP\tREF\tALT\tFILTER\tFather_GENO\tMother_GENO\tChild_GENO\tCODE\n")
    lmendel.write("CHR\tPOS\tSNP\tREF\tALT\tFILTER\tN\n")
    imendel.write("FID\tIID\tN\n")
    HetHaploid.write("FID\tIID\tCHR\tPOS\tSNP\tREF\tALT\tGENO\n")
    nMendelError = nHetHaploid = nMEtag = 0
    nchrY = nchrM = 0
    for line in input:
        if line.startswith("##"):
            if zeroout:
                output.write(line)
            continue
        elif line.startswith("#"): 
            header = line.strip().split()
            nIndiv = len(header) - 9
            niMendelError = [0 for i in header]
            print("{} individuals read from [ {} ]".format(nIndiv, vcf))
            vcfIDs = header[9:]
            if zeroout:
                output.write(line)
        elif line.startswith('Y') or line.startswith('chrY'):
            nchrY += 1
        elif line.startswith('M') or line.startswith('chrM'):
            nchrM += 1
        else:
            data = line.strip().split()
            chr = data[0]
            pos = data[1]
            rsID = data[2]
            refAllele = data[3]
            altAllele = data[4]
            filter = data[6]
            nlMendelError = 0
            #tranAlleles = str.maketrans('01',refAllele + altAllele,)
            if len(data[4].split(",")) != 1: # remove MAPs
                continue
            zeroout_index = []
            nTrioMendelError = 0
            if chr == 'X':
                for i in range(9,len(data)):
                    geno = data[i].split(":")[0]
                    info = data[i][3:]
                    if len(geno) == 3:
                        if SEX[header[i]] == '1' and geno.split(geno[1])[0] != geno.split(geno[1])[1]:
                            HetHaploid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(FID[header[i]],header[i],chr,pos,rsID,refAllele,altAllele,data[i]))
                            data[i] = './.' + info
                            nHetHaploid += 1
            for j in trios:
                father = header.index(PID[j])
                mother = header.index(MID[j])
                kid = header.index(j)
                kid_sex = SEX[j]
                fatherGeno = data[father].split(":")[0]
                fatherAlleles = set(fatherGeno.split(fatherGeno[1]))
                fatherAlleles.discard(".")
                motherGeno = data[mother].split(":")[0]
                motherAlleles = set(motherGeno.split(motherGeno[1]))
                motherAlleles.discard(".")
                kidGeno = data[kid].split(":")[0]
                kidAlleles = set(kidGeno.split(kidGeno[1]))
                kidAlleles.discard(".")
                parentsAlleles = fatherAlleles.union(motherAlleles)
                if len(fatherAlleles) == 1 and len(motherAlleles) == 1 and len(parentsAlleles) == 1 and len(kidAlleles) == 2:
                    nMendelError += 1
                    nlMendelError += 1
                    nTrioMendelError += 1
                    mendel.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(FID[j], header[father], header[mother], j, chr, pos, rsID, refAllele, altAllele, filter, data[father], data[mother], data[kid], 1))
                    zeroout_index.extend([father, mother, kid])
                elif len(fatherAlleles) == 1 and len(motherAlleles) == 1 and len(parentsAlleles) == 1 and len(kidAlleles) == 1 and kidAlleles.isdisjoint(parentsAlleles):
                    nMendelError += 1
                    nlMendelError += 1
                    nTrioMendelError += 1
                    mendel.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(FID[j], header[father], header[mother], j, chr, pos, rsID, refAllele, altAllele, filter, data[father], data[mother], data[kid], 2))
                    zeroout_index.append(kid)
                elif len(fatherAlleles) == 1 and len(kidAlleles) == 1 and kidAlleles.isdisjoint(fatherAlleles):
                    if chr == 'X' and kid_sex == '1':
                        pass
                    else:
                        nMendelError += 1
                        nlMendelError += 1
                        nTrioMendelError += 1
                        mendel.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(FID[j], header[father], header[mother], j, chr, pos, rsID, refAllele, altAllele, filter, data[father], data[mother], data[kid], 3))
                        zeroout_index.extend([father, kid])
                elif len(motherAlleles) == 1 and len(kidAlleles) == 1 and kidAlleles.isdisjoint(motherAlleles):
                    nMendelError += 1
                    nlMendelError += 1
                    nTrioMendelError += 1
                    mendel.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(FID[j], header[father], header[mother], j, chr, pos, rsID, refAllele, altAllele, filter, data[father], data[mother], data[kid], 4))
                    zeroout_index.extend([mother, kid])
            for j in fatherKids:
                father = header.index(PID[j])
                kid = header.index(j)
                kid_sex = SEX[j]
                if chr == 'X' and kid_sex == '1':
                    pass
                else:
                    fatherGeno = data[father].split(":")[0]
                    fatherAlleles = set(fatherGeno.split(fatherGeno[1]))
                    fatherAlleles.discard(".")
                    kidGeno = data[kid].split(":")[0]
                    kidAlleles = set(kidGeno.split(kidGeno[1]))
                    kidAlleles.discard(".")
                    if len(fatherAlleles) == 1 and len(kidAlleles) == 1 and kidAlleles.isdisjoint(fatherAlleles):
                        nMendelError += 1
                        nlMendelError += 1
                        mendel.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(FID[j], header[father], 'NA', j, chr, pos, rsID, refAllele, altAllele, filter, data[father], 'NA', data[kid], 3))
                        zeroout_index.extend([father, kid])
            for j in motherKids:
                mother = header.index(MID[j])
                kid = header.index(j)
                motherGeno = data[mother].split(":")[0]
                motherAlleles = set(motherGeno.split(motherGeno[1]))
                motherAlleles.discard(".")
                kidGeno = data[kid].split(":")[0]
                kidAlleles = set(kidGeno.split(kidGeno[1]))
                kidAlleles.discard(".")
                if len(motherAlleles) == 1 and len(kidAlleles) == 1 and kidAlleles.isdisjoint(motherAlleles):
                    nMendelError += 1
                    nlMendelError += 1
                    mendel.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(FID[j], 'NA', header[mother], j, chr, pos, rsID, refAllele, altAllele, filter, 'NA', data[mother], data[kid], 4))
                    zeroout_index.extend([mother, kid])
            lmendel.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr, pos, rsID, refAllele, altAllele, filter, nlMendelError))
            for i in set(zeroout_index):
                niMendelError[i] += 1
            if zeroout:
                for i in set(zeroout_index):
                    data[i] = './.' + data[i][3:]
                if nTrioMendelError/len(trios) > me:
                    if data[6] == '.' or data[6] == 'PASS':
                        data[6] = 'HighMendelianError'
                    else:
                        data[6] = data[6] + ';' + 'HighMendelianError'
                    output.write("{}\n".format("\t".join(j for j in data)))
                    nMEtag += 1
                else:
                    output.write("{}\n".format("\t".join(j for j in data)))
    for i in range(9,len(data)):
        imendel.write("{}\t{}\t{}\n".format(FID[header[i]], header[i], niMendelError[i]))
    print("Found {} heterozygous haploid genotypes; set to missing".format(nHetHaploid))
    print("Writing list of heterozygous haploid genotypes to [ {} ]".format(prefix + '.hh'))
    print("Found {} Mendelian Errors".format(nMendelError))
    print("Writing list of Mendelian Errors to [ {} ]".format(prefix + '.mendel'))
    if zeroout:
        print("Creating a new vcf [ {} ] by zeroing out all Mendelian Errors".format(prefix + '.vcf'))
        print("Marked {} variants as 'HighMendelianError' in the FILTER column".format(nMEtag))
        if nchrY >0 or nchrM > 0:
            print("Warning: {} variants on ChrY and {} variants on chrM were NOT included in the [ {} ]".format(nchrY, nchrM, prefix + '.vcf'))
    print()
    print("Analysis finished: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()

vcfPedcheck(args.vcf, args.fam, args.prefix, args.zeroout, args.me)
