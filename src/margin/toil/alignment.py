import os
import pysam
from toil_lib import require


class AlignmentFormat:
    SAM, BAM = range(2)


class AlignmentStruct(object):
    def __init__(self, filestoreId, alignment_format):
        # type (str, AlignmentFormat)
        self.fid        = filestoreId
        self.aln_format = alignment_format

    def FileStoreID(self):
        return self.fid

    def FID(self):
        return self.fid

    def AlignmentFormat(self):
        return self.aln_format

    def IsSAM(self):
        return (self.aln_format == AlignmentStruct.SAM)

    def IsBAM(self):
        return (self.aln_format == AlignmentFormat.BAM)


def splitLargeAlignment(parent_job, config, input_sam_fid):
    """takes a large alignment (SAM/BAM) and makes a bunch of smaller ones from it
    """
    def makeSamfile(batch):
        require(len(batch) > 0, "[splitLargeAlignment]This batch is empty")
        parent_job.fileStore.logToMaster("[splitLargeAlignment]Packing up {} alignments".format(len(batch)))
        temp_sam  = parent_job.fileStore.getLocalTempFileName()
        small_sam = pysam.Samfile(temp_sam, 'wh', template=sam)
        for aln in batch:
            small_sam.write(aln)
        small_sam.close()
        # write to JobStore
        fid = parent_job.fileStore.writeGlobalFile(temp_sam)
        return fid

    large_sam = parent_job.fileStore.readGlobalFile(input_sam_fid)
    require(os.path.exists(large_sam), "[splitLargeAlignment]Did't download large alignment")
    sam            = pysam.Samfile(large_sam, 'r')  # the big alignment
    small_sam_fids = []                             # list of FileStoreIDs of smaller alignments
    batch_of_alns  = []                             # list of alignedSegments
    total_alns     = 0                              # total alignments in the orig. to keep track
    for alignment in sam:
        if len(batch_of_alns) < config["split_alignments_to_this_many"]:  # add it to the batch
            batch_of_alns.append(alignment)
            total_alns += 1
            continue
        else:
            # write the alignedSegments in this batch to a new Samfile
            small_sam_fids.append(makeSamfile(batch_of_alns))
            # start a new batch and add this one
            batch_of_alns = []
            batch_of_alns.append(alignment)
            total_alns += 1

    if batch_of_alns != []:
        small_sam_fids.append(makeSamfile(batch_of_alns))

    parent_job.fileStore.logToMaster("[splitLargeAlignment]Input alignment has {n} alignments in it"
                                     "split it into {l} smaller files".format(n=total_alns,
                                                                              l=len(small_sam_fids)))
    return small_sam_fids
