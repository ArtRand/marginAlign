import os
import uuid
import pysam
import numpy as np
import pandas as pd
from collections import namedtuple
from toil_lib import require
from toil_lib.programs import docker_call
from margin.utils import getFastaDictionary
from margin.toil.localFileManager import LocalFile, urlDownload


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


def splitLargeAlignment(parent_job, split_alignments_to_this_many, input_sam_fid):
    """takes a large alignment (SAM/BAM) and makes a bunch of smaller ones from it
    """
    def makeSamfile(batch):
        require(len(batch) > 0, "[splitLargeAlignment]This batch is empty")
        parent_job.fileStore.logToMaster("[splitLargeAlignment]Packing up {} alignments".format(len(batch)))
        temp_sam  = parent_job.fileStore.getLocalTempFileName()
        small_sam = pysam.Samfile(temp_sam, 'wh', template=sam)
        for aln in batch:
            # make a UID for each alignment so we can look them up uniquely later
            uid              = uuid.uuid4().hex
            uid_to_read[uid] = aln.query_name
            aln.query_name   = uid
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
    uid_to_read    = {}                             # keep a map of uids for each alignment to reads
    total_alns     = 0                              # total alignments in the orig. to keep track

    for alignment in sam:
        if len(batch_of_alns) < split_alignments_to_this_many:  # add it to the batch
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
    return small_sam_fids, uid_to_read


def makeRangedAlignmentJobFunction(job, contig_label, full_alignment_fid, ranges,
                                   samtools_image="quay.io/ucsc_cgl/samtools"):
    uid               = uuid.uuid4().hex
    workdir           = job.fileStore.getLocalTempDir()
    ranged_alignments = []
    alignment         = LocalFile(workdir=workdir, filename="{r}{u}.bam".format(r=contig_label, u=uid))

    job.fileStore.readGlobalFile(full_alignment_fid, userPath=alignment.fullpathGetter())
    make_bai(job, alignment, workdir)

    for start, end in ranges:
        uid  = uuid.uuid4().hex
        aln  = LocalFile(workdir=workdir, filename="tmp{}.sam".format(uid))
        fH   = open(aln.fullpathGetter(), "w")
        args = ["view",
                "-h",
                "/data/{}".format(alignment.filenameGetter()),
                "{chr}:{start}-{end}".format(chr=contig_label, start=start, end=end - 1)]  # -1 bc. samtools
        docker_call(job=job, tool=samtools_image, parameters=args, work_dir=workdir, outfile=fH)
        fH.close()
        if check_for_empty(aln.fullpathGetter()):
            continue
        fid = job.fileStore.writeGlobalFile(aln.fullpathGetter())
        ranged_alignments.append(AlignmentShard(start=start, end=end, FileStoreID=fid))

    return ranged_alignments


def check_for_empty(alignment_file_path):
    "Returns True if empty"
    sam = pysam.Samfile(alignment_file_path, "r")
    for ar in sam:
        return False
    return True


def make_bai(parent_job, alignment_file, workdir, samtools_image="quay.io/ucsc_cgl/samtools"):
    require(os.path.exists(alignment_file.fullpathGetter()), "[make_bai]BAM file missing")
    make_bai_args = ["index", "/data/{}".format(alignment_file.filenameGetter())]
    docker_call(job=parent_job, tool=samtools_image, parameters=make_bai_args, work_dir=(workdir + "/"))


def get_ranges(reference_length, split_len):
    s      = 0
    e      = split_len
    step   = split_len
    ranges = []
    batch  = 10
    while e < reference_length:
        ranges.append((s, e))
        s += step
        e += step
    ranges.append((s, reference_length))
    return [ranges[i : i + batch] for i in xrange(0, len(ranges), batch)]


def shardAlignmentByRegionJobFunction(job, reference_fid, input_alignment_fid, chromosome_split_len,
                                      samtools_image="quay.io/ucsc_cgl/samtools"):
    def contig_generator():
        _args  = ["idxstats", "/data/{}".format(full_alignment.filenameGetter())]
        _stats = LocalFile(workdir=workdir, filename="{}.tmp".format(uuid.uuid4().hex))
        _fH    = open(_stats.fullpathGetter(), "w")
        docker_call(job=job,
                    tool=samtools_image,
                    parameters=_args,
                    work_dir=(workdir + "/"),
                    outfile=_fH)
        _fH.close()
        _table = pd.read_table(_stats.fullpathGetter(),
                               usecols=(0, 2),
                               names=["contig", "n_reads"],
                               dtype={"contig": np.str, "n_reads": np.int64})
        for contig in _table.loc[_table["n_reads"] > 0]["contig"]:
            yield contig

    uid            = uuid.uuid4().hex
    workdir        = job.fileStore.getLocalTempDir()
    reference_hash = getFastaDictionary(job.fileStore.readGlobalFile(reference_fid))
    full_alignment = LocalFile(workdir=workdir, filename="full{}.bam".format(uid))
    accumulator    = []

    job.fileStore.readGlobalFile(input_alignment_fid, userPath=full_alignment.fullpathGetter())
    make_bai(job, full_alignment, workdir)

    for contig in contig_generator():
        job.fileStore.logToMaster("[shardAlignmentByRegionJobFunction] %s contains alignments, sharding" % contig)
        contig_ranges  = get_ranges(len(reference_hash[contig]), chromosome_split_len)
        accumulator.extend([job.addChildJobFn(makeRangedAlignmentJobFunction,
                                              contig, input_alignment_fid, batch).rv()
                            for batch in contig_ranges])

    return accumulator


def downloadSplitAlignmentByRegionJobFunction(job, alignment_url, reference_fid, chromosome_split_len):
    """DEPRECIATED?
    """
    full_alignment = LocalFile(workdir=job.fileStore.getLocalTempDir(),
                               filename="{}.bam".format(uuid.uuid4().hex))
    urlDownload(parent_job=job, source_url=alignment_url, destination=full_alignment)
    require(os.path.exists(full_alignment.fullpathGetter()),
            "[downloadSplitAlignmentByRegionJobFunction]Didn't download alignment {}".format(alignment_url))
    alignment_fid = job.fileStore.writeGlobalFile(full_alignment.fullpathGetter())
    sharded_alignments = job.addChildJobFn(shardAlignmentByRegionJobFunction,
                                           reference_fid,
                                           alignment_fid,
                                           chromosome_split_len).rv()
    return alignment_fid, sharded_alignments


def shardAlignmentByRegion(parent_job, config, input_alignment_fid,
                           samtools_image="quay.io/ucsc_cgl/samtools"):
    """DEPRECIATED
    """
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
