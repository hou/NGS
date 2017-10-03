# Author: Liping Hou
# Email:  houliping@gmail.com
# Sept 18th, 2017

import argparse
import datetime 
import sys
import random
import numpy as np

parser = argparse.ArgumentParser(description="Perform permutation based gene-set enriment analysis")
parser.add_argument('genes', help="genes to be tested")
parser.add_argument('background', help='background gene list')
parser.add_argument('genesets', help="gene sets")
parser.add_argument('output', help='name of the result file')
parser.add_argument('--permutation', help='number of permutations (default: 10,000)', metavar="N", type=int, default=10000)
parser.add_argument('--match', help='match by the column "var" in the provided file',nargs=2, metavar=("file", "var"))
args = parser.parse_args()

def enrich(genes, background, genesets, outfile, perm, match):
    print()
    print("Analysis started: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()
    
    geneFile = open(genes)
    genesetFile = open(genesets)
    backgroundFile = open(background)
    output = open(outfile, 'w')

    nGeneSets = nGenes = nBackground = 0

    print("Reading genes from [ {} ]".format(genes))
    Genes = []
    for line in geneFile:
        Genes.append(line.strip())
        nGenes += 1
        
    TestGenes = set(Genes)
    print("{} lines (including {} unique genes) to be included from [ {} ]".format(nGenes, len(TestGenes), genes))

    print("Reading background genes from [ {} ]".format(background))
    BG_Genes = []
    for line in backgroundFile:
        BG_Genes.append(line.strip())
        nBackground += 1
        
    BackgroundGenes = set(BG_Genes)
    print("{} lines (including {} unique genes) to be included from [ {} ]".format(nBackground, len(BackgroundGenes), background))

    # Only genes found in background will be used for enrichment analysis
    TestGenesInBackground = TestGenes.intersection(BackgroundGenes)
    print("{} genes from [ {} ] overlap with background genes".format(len(TestGenesInBackground),genes))

    if match is not None:
        print("Reading '{}' from [ {} ]".format(match[1],match[0]))
        matchFile = open(match[0])
        matchHead = next(matchFile).strip().split()
        try:
            matchIndex = matchHead.index(match[1])
        except ValueError:
            print("Error: column '{}' does not exist in [ {} ]".format(match[1],match[0]))
            sys.exit()
        matchVal = []
        for line in matchFile:
            tem = line.strip().split()
            if tem[0] in BackgroundGenes:
                matchVal.append(tem[matchIndex])

        matchVal = [float(i) for i in matchVal]
        decileBreaks = np.percentile(np.array(matchVal),np.arange(0,101,10)).tolist()

        def decile(x, breaks):
            n = len(breaks) - 1
            i = 0
            while i < n:
                if x > breaks[i] and x <= breaks[i+1]:
                    return 'decile' + str(i+1)
                i += 1
            if x == breaks[0]: 
                return 'decile1'

        matchFile.seek(0)
        next(matchFile)
        DECILE = {}
        for i in range(1,11):
            DECILE['decile'+str(i)] = set()

        BackGroundGenesInMatch = set()
        for line in matchFile:
            tem = line.strip().split()
            if tem[0] in BackgroundGenes:
                grp = decile(float(tem[matchIndex]),decileBreaks)
                DECILE[grp].add(tem[0])
                BackGroundGenesInMatch.add(tem[0])

        if len(BackgroundGenes) - len(BackGroundGenesInMatch) >0:
            print("Warning: can not find '{}' for {} background genes".format(match[1], len(BackgroundGenes) - len(BackGroundGenesInMatch)))

        nDecile = [len(DECILE['decile'+str(i)]) for i in range(1,11)]
        print("Number of background genes in each decile: {}".format(" ".join(str(i) for i in nDecile)))

    output.write("Gene_set_name\tGene_set_size\tObserved_in_background\tObserved_in_tested\tExpected_in_tested\tPermutation\tPvalue\n")
    if match is not None:
        nTestGenesDecile = [] 
        TestGenesDecile = {}
        for i in range(1,11):
            TestGenesDecile['decile'+str(i)] = TestGenesInBackground.intersection(DECILE['decile'+str(i)])
            nTestGenesDecile.append(len(TestGenesDecile['decile'+str(i)]))
        print("Number of testing genes in each decile: {}".format(" ".join(str(i) for i in nTestGenesDecile)))

    for line in genesetFile:
        res = line.strip().split(maxsplit=1)
        setName = res[0]
        GeneNames = res[1].split(",")
        GeneNames = set([i.strip(" ") for i in GeneNames])
        n = len(GeneNames)
        counts = []
        if match is not None:
            nInBackground = len(GeneNames.intersection(BackGroundGenesInMatch))
            nObs = sum([len(GeneNames.intersection(TestGenesDecile['decile'+str(i)])) for i in range(1,11)])
            nInDecile = [len(GeneNames.intersection(DECILE['decile'+str(j)])) for j in range(1,11)]
            print("{} - number of genes in each decile: {}".format(setName," ".join(str(i) for i in nInDecile)))
            nExp = sum([nInDecile[i]/nDecile[i]*nTestGenesDecile[i] for i in range(0,10)])
            for j in range(perm):
                countsDecile = [len(set(random.sample(DECILE['decile'+str(k)], nTestGenesDecile[k-1])).intersection(GeneNames)) for k in range(1,11)]
                counts.append(sum(countsDecile))
        else:
            nInBackground = len(GeneNames.intersection(BackgroundGenes))
            nObs = len(TestGenesInBackground.intersection(GeneNames))
            nExp = nInBackground/len(BackgroundGenes)*len(TestGenesInBackground)
            for i in range(perm):
                randomGenes = random.sample(BackgroundGenes, len(TestGenesInBackground))
                counts.append(len(set(randomGenes).intersection(GeneNames)))
        GTn = len([i for i in counts if i > nObs])
        pval = (1+GTn)/(perm + 1)
        output.write("{}\t{}\t{}\t{}\t{:.2f}\t{}\t{:.2E}\n".format(setName, n, nInBackground, nObs, nExp, perm, pval))
        
    geneFile.close()
    genesetFile.close()
    backgroundFile.close()
    if match is not None:
        matchFile.close()
    output.close()

    print("Output file was written to: [ {} ]".format(outfile))
    print()
    print("Analysis finished: {}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print()

enrich(args.genes, args.background, args.genesets, args.output, args.permutation, args.match)
