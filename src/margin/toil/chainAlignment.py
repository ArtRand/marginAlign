"""Functions to chain an alignment
"""
import pysam
from margin.utils import \
    getFirstNonClippedPositionInRead, \
    getLastNonClippedPositionInRead, \
    getFastaDictionary, \
    samIterator

from sonLib.bioio import reverseComplement, fastqRead


def mergeChainedAlignedSegments(chainedAlignedSegments, refSequence, readSequence):
    """Makes a single alignment for the given chained reads. Will soft soft clip
    the unaligned prefix and suffix of the readSequence.

    From doc on building pysam line
    a = pysam.AlignedSegment()
    a.query_name = "read_28833_29006_6945"
    a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 32
    a.mapq = 20
    a.cigar = ( (0,10), (2,1), (0,25) )
    a.mrnm = 0
    a.mpos=199
    a.isize=167
    a.qual="<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
    a.tags = ( ("NM", 1),
               ("RG", "L1") )
    """
    cAR = pysam.AlignedSegment()
    aR = chainedAlignedSegments[0]
    cAR.query_name = aR.query_name

    # Parameters we don't and therefore set properly
    # cAR.flag = aR.flag
    # cAR.mapq = aR.mapq
    # cAR.mrnm = 0
    # cAR.mpos=0
    # cAR.isize=0
    # cAR.qual = "<" * len(readSequence)
    # cAR.tags = aR.tags 
    cAR.next_reference_id = -1
    cAR.reference_start = aR.reference_start  # Reference start
    cAR.is_reverse = aR.is_reverse
    cAR.query_sequence = reverseComplement(readSequence) if cAR.is_reverse else readSequence
    cAR.reference_id = aR.reference_id
    cigarList = []
    pPos = aR.reference_start
    # Iterate from the other end of the sequence if reversed
    pQPos = -(len(readSequence) - 1) if cAR.is_reverse else 0

    for aR in chainedAlignedSegments:
        assert cAR.is_reverse == aR.is_reverse
        # Add a deletion representing the preceding unaligned reference positions
        assert aR.reference_start >= pPos
        if aR.reference_start > pPos:
            cigarList.append((2, aR.reference_start - pPos))
            pPos = aR.reference_start

        # Add an insertion representing the preceding unaligned read positions
        # make it a soft clip if it is the first chained alignment
        qPos = getFirstNonClippedPositionInRead(aR, readSequence)
        assert qPos >= pQPos
        if qPos > pQPos:
            cigarList.append((4 if aR == chainedAlignedSegments[0] else 1, qPos - pQPos))
            pQPos = qPos

        # Add the operations of the cigar, filtering hard and soft clipping
        for op, length in aR.cigar:
            assert op in (0, 1, 2, 4, 5)
            if op in (0, 1, 2):
                cigarList.append((op, length))
            if op in (0, 2):  # Is match or deletion
                pPos += length
            if op in (0, 1):  # Is match or insertion
                pQPos += length

    assert pPos <= len(refSequence)

    # Set reference end coordinate (which is exclusive)
    # cAR.reference_end = pPos #We don't do this because it is set by cigar string

    # Now add any trailing, necessary soft clipping
    if cAR.is_reverse:
        assert pQPos <= 1
        if pQPos < 1:
            cigarList.append((4, -pQPos + 1))
    else:
        assert pQPos <= len(readSequence)
        if pQPos < len(readSequence):
            cigarList.append((4, len(readSequence) - pQPos))

    cAR.cigar = tuple(cigarList)

    # Check ops
    for op, length in cAR.cigar:  # We should have no hard clipped ops
        assert op in (0, 1, 2, 4)

    # Reference sequence check coordinates
    assert sum([length for op, length in cigarList if op in (0, 2)]) == cAR.reference_end - cAR.reference_start
    assert cAR.reference_start >= 0 and cAR.reference_start < len(refSequence)
    assert cAR.reference_end >= 0 and cAR.reference_end <= len(refSequence)

    # Read sequence check coordinates
    assert cAR.query_alignment_start >= 0 and cAR.query_alignment_start < len(readSequence)
    assert cAR.query_alignment_end >= 0 and cAR.query_alignment_end <= len(readSequence)
    assert cAR.query_alignment_start + sum([length for op, length in cigarList if op in (0, 1)]) == cAR.query_alignment_end

    return cAR


def chainFn(alignedSegments, refSeq, readSeq,
            scoreFn=lambda alignedSegment, refSeq, readSeq : sum([length for op, length in alignedSegment.cigar if op == 0]),
            maxGap=200):
    # Score function is number of aligned pairs
    """Gets the highest scoring chain of alignments on either the forward or reverse
    strand. Score is (by default) number of aligned positions.
    """
    def getStartAndEndCoordinates(alignedSegment):
        """Gets the start and end coordinates in both the reference and query, using coordinates
        relative to the original read and reference equence
        """
        return alignedSegment.reference_start, getFirstNonClippedPositionInRead(alignedSegment, readSeq), \
            alignedSegment.reference_end - 1, getLastNonClippedPositionInRead(alignedSegment, readSeq)

    alignedSegmentToScores = dict([(aR, scoreFn(aR, refSeq, readSeq)) for aR in alignedSegments])
    alignedSegmentToCoordinates = dict([(aR, getStartAndEndCoordinates(aR)) for aR in alignedSegments])
    alignedSegmentPointers = {}

    # Currently uses sloppy quadratic algorithm to find highest chain
    alignedSegments = sorted(alignedSegments, key=lambda aR : alignedSegmentToCoordinates[aR][0])
    # Sort by reference coordinate
    for i in xrange(len(alignedSegments)):
        aR = alignedSegments[i]
        rStart, qStart, rEnd, qEnd = alignedSegmentToCoordinates[aR]
        score = alignedSegmentToScores[aR]
        for j in xrange(i):  # Look at earlier alignments in list
            aR2 = alignedSegments[j]
            rStart2, qStart2, rEnd2, qEnd2 = alignedSegmentToCoordinates[aR2]
            assert rStart2 <= rStart
            # conditions for a chain
            if (rStart > rEnd2 and
               qStart > qEnd2 and
               aR.is_reverse == aR2.is_reverse and
               rStart - rEnd2 + qStart - qEnd2 <= maxGap and
               score + alignedSegmentToScores[aR2] > alignedSegmentToScores[aR]):
                alignedSegmentToScores[aR] = score + alignedSegmentToScores[aR2]
                alignedSegmentPointers[aR] = aR2

    # Now find highest scoring alignment 
    aR = max(alignedSegments, key=lambda aR : alignedSegmentToScores[aR])

    # Construct chain of alignedSegments
    chain = [aR]
    while aR in alignedSegmentPointers:
        aR = alignedSegmentPointers[aR]
        chain.append(aR)
    chain.reverse()

    return chain


def chainSamFile(parent_job, samFile, outputSamFile, readFastqFile, referenceFastaFile,
                 chainFn=chainFn):
    """Chains together the reads in the SAM file so that each read is covered by a
    single maximal alignment
    """
    sam = pysam.Samfile(samFile, "r")
    refSequences = getFastaDictionary(referenceFastaFile)  # Hash of names to sequences

    alignmentsHash = {}
    for aR in samIterator(sam):  # Iterate on the sam lines and put into buckets by read
        # This should be improved, because the whole sam file is being stored in memory
        if aR.query_name not in alignmentsHash:
            alignmentsHash[aR.query_name] = {}
        if aR.reference_id not in alignmentsHash[aR.query_name]:
            alignmentsHash[aR.query_name][aR.reference_id] = []
        alignmentsHash[aR.query_name][aR.reference_id].append(aR)

    # Now write out the sam file
    outputSam = pysam.Samfile(outputSamFile, "wb", template=sam)

    # Chain together the reads
    chainedAlignedSegments = []
    for readName, readSeq, qualValues in fastqRead(readFastqFile):
        readName = readName.split()[0]  # Remove any white space from the name
        if readName in alignmentsHash:
            for refID in alignmentsHash[readName].keys():
                alignedSegments = alignmentsHash[readName][refID]
                try:
                    refSeq = refSequences[sam.getrname(refID)]
                    # XXX could make this a child function instead and allow the chaining to be done in ||
                    chainedAlignedSegments.append(mergeChainedAlignedSegments(chainFn(alignedSegments,
                                                  refSeq, readSeq), refSeq, readSeq))
                except KeyError:
                    parent_job.fileStore.logToMaster("[chainSamFile]Missing reference sequence {}"
                                                     "".format(sam.getrname(refID)))
                    continue
            alignmentsHash.pop(readName)

    if len(alignmentsHash) != 0:
        parent_job.fileStore.logToMaster("[chainSamFile]WARNING not all reads were re-chained")

    # Sort chained alignments by reference and reference coordinates
    chainedAlignedSegments.sort(key=lambda aR : (sam.getrname(aR.reference_id),
                                                 aR.reference_start, aR.reference_end))

    for cAR in chainedAlignedSegments:
        outputSam.write(cAR)

    sam.close()
    outputSam.close()
