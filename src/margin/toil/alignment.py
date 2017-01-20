import os
import uuid
import pysam
from collections import namedtuple
from toil_lib import require
from toil_lib.programs import docker_call
from margin.utils import getFastaDictionary
from margin.toil.localFileManager import LocalFile


class AlignmentFormat:
    SAM, BAM = range(2)

AlignmentShard = namedtuple("AlignmentShard", ["start", "end", "FileStoreID"])


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
        return AlignmentShard(start=None, end=None, FileStoreID=fid)

    large_sam = parent_job.fileStore.readGlobalFile(input_sam_fid)
    require(os.path.exists(large_sam), "[splitLargeAlignment]Did't download large alignment")
    sam            = pysam.Samfile(large_sam, 'r')  # the big alignment
    small_sam_fids = []                             # list of FileStoreIDs of smaller alignments
    batch_of_alns  = []                             # batch of alignedSegments
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

    parent_job.fileStore.logToMaster("[splitLargeAlignment]Input alignment has {N} alignments in it"
                                     "split it into {n} smaller files".format(N=total_alns,
                                                                              n=len(small_sam_fids)))
    return small_sam_fids


def shardAlignmentByRegion(parent_job, config, input_alignment_fid,
                           samtools_image="quay.io/ucsc_cgl/samtools"):
    def get_ranges(reference_length):
        s      = 0
        e      = config["split_chromosome_this_length"]
        step   = config["split_chromosome_this_length"]
        ranges = []
        while e < reference_length:
            ranges.append((s, e))
            s += step
            e += step
        ranges.append((s, reference_length))
        return ranges

    def check_for_empty(alignment_file_path):
        "Returns True if empty"
        sam = pysam.Samfile(alignment_file_path, "r")
        for ar in sam:
            return False
        return True

    def make_bai(alignment_file):
        require(os.path.exists(alignment_file.fullpathGetter()), "[make_bai]BAM file missing")
        make_bai_args = ["index", "/data/{}".format(alignment_file.filenameGetter())]
        docker_call(job=parent_job, tool=samtools_image, parameters=make_bai_args, work_dir=(workdir + "/"))

    def break_alignment_by_region(alignment_file, contig_label):
        reference_length  = len(reference_hash[contig_label])
        ranged_alignments = []
        make_bai(alignment_file)
        # loop over the ranges and make sub-alignments
        for start, end in get_ranges(reference_length):
            require(os.path.exists(alignment_file.fullpathGetter()), "[break_alignment_by_region]DIdn't find SAM")
            aln  = parent_job.fileStore.getLocalTempFile()
            fH   = open(aln, "w")
            args = ["view",
                    "-h",
                    "/data/{}".format(alignment_file.filenameGetter()),
                    "{chr}:{start}-{end}".format(chr=contig_label, start=start, end=end - 1)]  # -1 bc. samtools
            docker_call(job=parent_job, tool=samtools_image, parameters=args, work_dir=(workdir + "/"), outfile=fH)
            fH.close()
            if check_for_empty(aln):
                continue
            fid = parent_job.fileStore.writeGlobalFile(aln)
            ranged_alignments.append(AlignmentShard(start=start, end=end, FileStoreID=fid))

        return ranged_alignments

    # an unfortunate reality is that we need to check every contig in the input reference for an alignment
    # even though in practice we'll have an alignment sorted to contain only one contig
    uid            = uuid.uuid4().hex
    workdir        = parent_job.fileStore.getLocalTempDir()
    reference_hash = getFastaDictionary(parent_job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
    full_alignment = LocalFile(workdir=workdir, filename="full{}.bam".format(uid))
    accumulator    = []  # will contain [(start, end, fid)...] for each alignment shard
    parent_job.fileStore.readGlobalFile(input_alignment_fid, userPath=full_alignment.fullpathGetter())
    make_bai(full_alignment)  # make the BAI for the fill alignment

    # loop over the contigs in our reference 
    for reference in reference_hash:
        contig_alignment = LocalFile(workdir=workdir,
                                     filename="{contig}{uid}.bam".format(contig=reference, uid=uid))
        samtools_args    = ["view",
                            "-hb",
                            "/data/{}".format(full_alignment.filenameGetter()),
                            "{}".format(reference)]
        _handle          = open(contig_alignment.fullpathGetter(), "w")
        docker_call(job=parent_job,
                    tool=samtools_image,
                    parameters=samtools_args,
                    work_dir=(workdir + "/"),
                    outfile=_handle)
        _handle.close()
        if check_for_empty(contig_alignment.fullpathGetter()):
            continue
        accumulator.extend(break_alignment_by_region(contig_alignment, reference))
 
    return accumulator
