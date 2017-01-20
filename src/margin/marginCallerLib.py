import vcf
import datetime
import math

from itertools import product
from margin.utils import *
from margin.toil.hmm import Hmm

BASES = "ACGT"


def getProb(subMatrix, start, end):
    """Get the substitution probability.
    """
    return subMatrix[(start, end)]


def getEvoProb(evoMatrix, start, end):
    if start in BASES:
        return evoMatrix[(start, end)]
    else:
        return 0.25


def calcBasePosteriorProbs(baseObservations, refBase,
                           evolutionarySubstitionMatrix, errorSubstutionMatrix):
    """Function that does the column probability calculation.
    """
    logBaseProbs = map(lambda missingBase : \
            math.log(getEvoProb(evolutionarySubstitionMatrix, refBase.upper(), missingBase)) + 
            reduce(lambda x, y : x + y, map(lambda observedBase : \
                        math.log(getProb(errorSubstutionMatrix, missingBase, 
                                         observedBase))*baseObservations[observedBase], BASES)), BASES)
    totalLogProb = reduce(lambda x, y : x + math.log(1 + math.exp(y-x)), logBaseProbs)
    return dict(zip(BASES, map(lambda logProb : math.exp(logProb - totalLogProb), logBaseProbs)))

def loadHmmSubstitutionMatrix(hmmFile):
    """Load the substitution matrix from an HMM file
    """
    hmm = Hmm.loadHmm(hmmFile)
    m = hmm.emissions[:len(BASES)**2]
    m = map(lambda i : m[i] / sum(m[4*(i/4):4*(1 + i/4)]), range(len(m))) #Normalise m
    return dict(zip(product(BASES, BASES), m))

def getNullSubstitutionMatrix():
    """Null matrix that does nothing
    """
    return dict(zip(product(BASES, BASES), [1.0]*len(BASES)**2))

def vcfRead(vcfFile):
    vcfCalls = set()
    for x in vcf.Reader(open(vcfFile, 'r')):
        for alt in x.ALT:
            vcfCalls.add((x.CHROM, int(x.POS), str(alt).upper()))
    return vcfCalls

def vcfWriteHeader(outfile, referenceFastaFile):
    vcfFile = open(outfile, "w")
    vcfFile.write("##fileformat=VCFv4.2\n")
    vcfFile.write("##fileDate=" + str(datetime.datetime.now().date()).replace("-", "") + "\n")
    vcfFile.write("##source=marginCaller\n")
    vcfFile.write("##reference=" + referenceFastaFile + "\n")
    vcfFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
    vcfFile.close()
    return

def vcfWrite(referenceFastaFile, refSequences, variantCalls, outputVcfFile, write_header=True):
    #Organise the variant calls
    variantCallsHash = dict(map(lambda x : (x, {}), refSequences))
    for refSeqName, refPosition, base, posteriorProb in variantCalls:
        if not refPosition in variantCallsHash[refSeqName]:
            variantCallsHash[refSeqName][refPosition] = []
        variantCallsHash[refSeqName][refPosition].append((base, posteriorProb))

    if write_header:
        vcfWriteHeader(outputVcfFile, referenceFastaFile)
    vcfFile = open(outputVcfFile, "w")
    # iterate through records now
    for refSeqName in refSequences:
        for refPosition, refBase in enumerate(refSequences[refSeqName]):
            chrom = refSeqName
            pos = str(refPosition + 1) # 1-based indexing in reference
            id = "."
            ref = refBase
            qual = "."
            filter = "PASS"
            fmt = ""
            alt = "."
            info = "."
            
            if refPosition in variantCallsHash[refSeqName].keys():
                variant = []
                posteriorProb = []
                for variantCall in variantCallsHash[refSeqName][refPosition]:
                    variant.append(variantCall[0])
                    posteriorProb.append(str(variantCall[1]))
                alt = ",".join(variant)
                info = ",".join(posteriorProb)
                # _Record(chrom, pos, ID, ref, alt, qual, filt, info, fmt, self._sample_indexes)
                Record = refSeqName + "\t" + pos + "\t" + id + "\t" + ref + "\t" + \
                alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + fmt
                vcfFile.write(Record)
                vcfFile.write("\n")
    vcfFile.close()
