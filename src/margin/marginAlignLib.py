#import pysam
#import sys
#import os
#from jobTree.src.bioio import reverseComplement, fastaRead, system, fastaWrite, \
#cigarRead, logger, nameValue
from margin.utils import *
#from cPecan import cPecanEm
#from cPecan.cPecanEm import Hmm, SYMBOL_NUMBER
import numpy as np

"""
def learnModelFromSamFileTargetFn(target, samFile, readFastqFile, 
                                  referenceFastaFile, options):
    //Does expectation maximisation on sam file to learn the hmm for the sam file.
    
    #Get cigars and reads fasta file
    cigars = os.path.join(target.getGlobalTempDir(), "temp.cigar")
    fHCigars = open(cigars, 'w')
    reads = os.path.join(target.getGlobalTempDir(), "temp.fa")
    fHReads = open(reads, 'w')
    sam = pysam.Samfile(samFile, "r" )
    for aR, counter in zip(sam, xrange(sys.maxint)): #Iterate on the sam lines realigning them in parallel            
        aR.query_name = aR.query_name + "_%s" % counter
        fHCigars.write(getExonerateCigarFormatString(aR, sam) + "\n")
        fastaWrite(fHReads, aR.query_name, aR.seq)
    fHCigars.close(); fHReads.close()
    
    unnormalisedOutputModel = os.path.join(target.getGlobalTempDir(), 
                                           "unnormalisedOutputModel.hmm")
    target.addChildTargetFn(cPecanEm.expectationMaximisationTrials, 
                            args=(" ".join([reads, referenceFastaFile ]), cigars, 
                                  unnormalisedOutputModel, options))
    
    #Now set up normalisation
    target.setFollowOnTargetFn(learnModelFromSamFileTargetFn2, 
                               args=(unnormalisedOutputModel, options))

def learnModelFromSamFileTargetFn2(target, unnormalisedOutputModel, options):
    hmm = Hmm.loadHmm(unnormalisedOutputModel)
    setHmmIndelEmissionsToBeFlat(hmm)
    #Normalise background emission frequencies, if requested to GC% given
    normaliseHmmByReferenceGCContent(hmm, 0.5)
    hmm.write(options.outputModel)
    
    """
#toMatrix = lambda e : map(lambda i : e[SYMBOL_NUMBER*i:SYMBOL_NUMBER*(i+1)], 
#                          xrange(SYMBOL_NUMBER))
#fromMatrix = lambda e : reduce(lambda x, y : list(x) + list(y), e)
    
def normaliseHmmByReferenceGCContent(hmm, gcContent):
    """Normalise background emission frequencies to GC% given
    """
    for state in range(hmm.stateNumber):
        if state not in (2, 4): #Don't normalise GC content of insert states 
            #(as they don't have any ref bases!)
            n = toMatrix(hmm.emissions[(SYMBOL_NUMBER**2) * 
                                       state:(SYMBOL_NUMBER**2) * (state+1)])
            hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = \
            fromMatrix(map(lambda i : map(lambda j : (n[i][j]/sum(n[i])) * 
            (gcContent/2.0 if i in [1, 2] else (1.0-gcContent)/2.0), range(SYMBOL_NUMBER)), 
                           range(SYMBOL_NUMBER))) #Normalise

def setHmmIndelEmissionsToBeFlat(hmm):
    """Set indel emissions to all be flat
    """
    for state in range(1, hmm.stateNumber):
        hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = \
        [1.0/(SYMBOL_NUMBER**2)]*SYMBOL_NUMBER**2  

def modifyHmmEmissionsByExpectedVariationRate(hmm, substitutionRate):
    #Normalise background emission frequencies, if requested to GC% given
    n = toMatrix(map(lambda i : (1.0-substitutionRate) if i % SYMBOL_NUMBER == \
                     i / SYMBOL_NUMBER else substitutionRate/(SYMBOL_NUMBER-1), 
                     xrange(SYMBOL_NUMBER**2)))
    hmm.emissions[:SYMBOL_NUMBER**2] = fromMatrix(np.dot(toMatrix(hmm.emissions[:SYMBOL_NUMBER**2]), n))

"""
def realignSamFileTargetFn(target, samFile, outputSamFile, readFastqFile, 
                           referenceFastaFile, options, chainFn=chainFn):
    //Chains and then realigns the resulting global alignments, using jobTree to 
    do it in parallel on a cluster.
    Optionally runs expectation maximisation.
    
    #Optionally chain the sam file
    if not options.noChain:
        target.logToMaster("Going to chain sam file: %s" % samFile)
        tempSamFile = os.path.join(target.getGlobalTempDir(), "temp.sam")
        chainSamFile(samFile, tempSamFile, readFastqFile, referenceFastaFile, chainFn)
        samFile = tempSamFile
    
    #If we do expectation maximisation we split here:
    if options.em:
        target.logToMaster("Going to run EM training with sam file: %s, read fastq file: %s, \
        reference fasta file: %s" \
                           % (samFile, readFastqFile, referenceFastaFile))
        target.addChildTargetFn(learnModelFromSamFileTargetFn, args=(samFile, 
                                    readFastqFile, referenceFastaFile, options))

    options.hmmFile = options.outputModel if options.em else options.inputModel #This
    #setups the hmm to be used the realignment function
    
    target.logToMaster("Going to realign sam file: %s to create output sam file: %s \
    with match gamma %s and gap gamma %s and model %s" % (samFile, outputSamFile, 
                options.gapGamma, options.matchGamma, options.hmmFile))
    
    target.setFollowOnTargetFn(paralleliseSamProcessingTargetFn, 
                               args=(samFile, 
                                     referenceFastaFile, outputSamFile, 
                                     realignCigarTargetFn, realignSamFile3TargetFn,
                                     options))
"""
"""
def realignCigarTargetFn(target, exonerateCigarStringFile, referenceSequenceName, 
                         referenceSequence, querySequenceFile, 
                         outputCigarFile, options):
    #Temporary files
    tempRefFile = os.path.join(target.getLocalTempDir(), "ref.fa")
    tempReadFile = os.path.join(target.getLocalTempDir(), "read.fa")
    
    #Write the temporary reference file.
    fastaWrite(tempRefFile, referenceSequenceName, referenceSequence) 
    
    #For each cigar string
    for exonerateCigarString, (querySequenceName, querySequence) in \
    zip(open(exonerateCigarStringFile, "r"), fastaRead(querySequenceFile)):
        fastaWrite(tempReadFile, querySequenceName, querySequence)
        #Call to cPecanRealign
        loadHmm = nameValue("loadHmm", options.hmmFile)
        system("echo \"%s\" | cPecanRealign %s %s --diagonalExpansion=10 \
        --splitMatrixBiggerThanThis=3000 %s --gapGamma=%s --matchGamma=%s >> %s" % \
               (exonerateCigarString[:-1], tempRefFile, tempReadFile, loadHmm, 
                options.gapGamma, options.matchGamma, outputCigarFile))
"""
"""
def realignSamFile3TargetFn(target, samFile, referenceFastaFile, 
                            outputSamFile, tempCigarFiles, options):
    #Setup input and output sam files
    sam = pysam.Samfile(samFile, "r" )
    
    #Replace the cigar lines with the realigned cigar lines
    outputSam = pysam.Samfile(outputSamFile, "wh", template=sam)
    def cigarIterator():
        #Iterates over all the cigars in the temp files.
        for tempCigarFile in tempCigarFiles:
            for pA in cigarRead(open(tempCigarFile)):
                yield pA 
        yield None #This is put in to cause an error if there is fewer 
        #cigars than pairwise alignments
    for aR, pA in zip(samIterator(sam), cigarIterator()): #Iterate on the sam lines 
        #realigning them in parallel
        
        #Replace alignment by converting exonerate ops to aligned read ops,
        #adding soft clipping/hard clipping of unaligned prefix and suffix of read
        ops = [ ]
        if len(aR.cigar) > 0 and aR.cigar[0][0] == 5:
            ##Add any hard clipped prefix
            ops.append(aR.cigar[0])
        if aR.query_alignment_start > 0:
            ops.append((4, aR.qstart))
        ops += map(lambda op : (op.type, op.length), pA.operationList)
        if aR.query_alignment_end < len(aR.query_sequence):
            ops.append((4, len(aR.query_sequence) - aR.query_alignment_end))
        if len(aR.cigar) > 1 and aR.cigar[-1][0] == 5: 
            ##Add any hard clipped suffix
            ops.append(aR.cigar[-1])
        
        #Checks the final operation list 
        ##Correct for the read
        assert sum(map(lambda (type, length) : length if type in (0,1,4) else 0, ops)) == \
        sum(map(lambda (type, length) : length if type in (0,1,4) else 0, aR.cigar))
        ##Correct for the reference
        assert sum(map(lambda (type, length) : length if type in (0, 2) else 0, ops)) == \
        aR.reference_end - aR.reference_start
        
        aR.cigar = tuple(ops)
        
        #Write out
        outputSam.write(aR)
    
    #Finish up
    sam.close()
    outputSam.close()
    
    target.logToMaster("Realigned sam file: %s to create output sam file: %s" % \
                       (samFile, outputSamFile))
"""
