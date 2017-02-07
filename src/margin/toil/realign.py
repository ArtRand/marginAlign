"""Functions to realign a SAM file with the base-level HMM
"""
import os
import pysam
import uuid
import cPickle

from toil_lib import require
from toil_lib.programs import docker_call

from sonLib.bioio import cigarRead
from margin.toil.hmm import downloadHmm
from margin.utils import getAlignedSegmentDictionary

from alignment import splitLargeAlignment
from localFileManager import LocalFile, deliverOutput
from shardAlignment import shardSamJobFunction


DOCKER_DIR = "/data/"


def setupLocalFiles(parent_job, global_config, hmm):
    uid             = uuid.uuid4().hex  # uid to be super sure we don't have any file collisions
    workdir         = parent_job.fileStore.getLocalTempDir()
    local_hmm       = LocalFile(workdir=workdir, filename="hmm.%s.txt" % uid)
    local_output    = LocalFile(workdir=workdir, filename="cPecan_out.%s.txt" % uid)
    local_input_obj = LocalFile(workdir=workdir, filename="cPecanInput.%s.pkl" % uid)
    hmm.write(local_hmm.fullpathGetter())
    return workdir, local_hmm, local_output, local_input_obj


def sortResultsByBatch(cPecan_result_fids):
    batch_sorted       = sorted(cPecan_result_fids, key=lambda tup: tup[1])
    sorted_cPecan_fids = [x[0] for x in batch_sorted]
    return [x for x in sorted_cPecan_fids if x is not None]


def realignSamFileJobFunction(job, config, input_samfile_fid, output_label):
    smaller_alns, uid_to_read = splitLargeAlignment(job,
                                                    config["split_alignments_to_this_many"],
                                                    input_samfile_fid)
    realigned_fids      = []
    hidden_markov_model = downloadHmm(job, config)

    for aln in smaller_alns:
        disk   = input_samfile_fid.size + config["reference_FileStoreID"].size
        memory = (6 * input_samfile_fid.size)
        realigned_fids.append(job.addChildJobFn(shardSamJobFunction, config, aln, hidden_markov_model,
                                                cPecanRealignJobFunction,
                                                rebuildSamJobFunction,
                                                batch_disk=disk,
                                                followOn_disk=(2 * config["reference_FileStoreID"].size),
                                                followOn_mem=(6 * aln.FileStoreID.size),
                                                disk=disk, memory=memory).rv())

    job.addFollowOnJobFn(combineRealignedSamfilesJobFunction, config, input_samfile_fid, realigned_fids,
                         uid_to_read, output_label)


def combineRealignedSamfilesJobFunction(job, config, input_samfile_fid, realigned_fids, uid_to_read, output_label):
    original_sam      = job.fileStore.readGlobalFile(input_samfile_fid)
    sam               = pysam.Samfile(original_sam, "r")
    filename          = "{sample}_{out_label}.bam".format(sample=config["sample_label"], out_label=output_label)
    output_sam        = LocalFile(workdir=job.fileStore.getLocalTempDir(), filename=filename)
    output_sam_handle = pysam.Samfile(output_sam.fullpathGetter(), "wb", template=sam)
    sam.close()

    for fid in realigned_fids:
        local_copy = job.fileStore.readGlobalFile(fid)
        samfile    = pysam.Samfile(local_copy, "rb")
        for alignment in samfile:
            read_label = uid_to_read[alignment.query_name]
            alignment.query_name = read_label
            output_sam_handle.write(alignment)
        samfile.close()
        job.fileStore.deleteGlobalFile(fid)

    output_sam_handle.close()
    deliverOutput(job, output_sam, config["output_dir"])


def cPecanRealignJobFunction(job, global_config, job_config, hmm, batch_number,
                             cPecan_image="quay.io/artrand/cpecanrealign"):
    # type: (toil.job.Job, dict<string, parameters>, dict<string, string>)
    """Runs Docker-ized cPecan HMM
    """
    workdir, local_hmm, local_output, local_input_obj = setupLocalFiles(job, global_config, hmm)
    if global_config["debug"]:
        job.fileStore.logToMaster("[cPecanRealignJobFunction]Batch {batch} using HMM from {model} "
                                  "and EM is {em}".format(batch=batch_number,
                                                          model=local_hmm.filenameGetter(),
                                                          em=global_config["EM"]))

    # pickle the job_config, that contains the reference sequence, the query sequences, and 
    # the pairwise alignments in exonerate format
    _handle = open(local_input_obj.fullpathGetter(), "w")
    cPickle.dump(job_config, _handle)
    _handle.close()

    # run cPecan in a container
    input_arg         = "--input={}".format(DOCKER_DIR + local_input_obj.filenameGetter())
    hmm_arg           = "--hmm_file={}".format(DOCKER_DIR + local_hmm.filenameGetter())
    gap_gamma_arg     = "--gap_gamma={}".format(global_config["gap_gamma"])
    match_gamma_arg   = "--match_gamma={}".format(global_config["match_gamma"])
    output_arg        = "--output_alignment_file={}".format(DOCKER_DIR + local_output.filenameGetter())
    cPecan_parameters = [input_arg, hmm_arg, gap_gamma_arg, match_gamma_arg, output_arg]
    try:
        docker_call(job=job,
                    tool=cPecan_image,
                    parameters=cPecan_parameters,
                    work_dir=(workdir + "/"))
        if os.path.exists(local_output.fullpathGetter()):
            return job.fileStore.writeGlobalFile(local_output.fullpathGetter())
        else:
            return None
    except:
        return None


def rebuildSamJobFunction(job, config, alignment_shard, cPecan_cigar_fileIds):
    # iterates over the files, downloads them, and iterates over the alignments 
    def cigar_iterator():
        for fid in sorted_cPecan_fids:
            local_copy = job.fileStore.readGlobalFile(fid)
            require(os.path.exists(local_copy),
                    "[cigar_iterator]ERROR couldn't find alignment {}".format(local_copy))
            for pA in cigarRead(open(local_copy, 'r')):
                yield pA
            job.fileStore.deleteLocalFile(fid)

    # get the chained alignment, we're going to replace the the alignments with the cPecan alignments
    alignment_fid  = alignment_shard.FileStoreID
    local_sam_path = job.fileStore.readGlobalFile(alignment_fid)

    try:
        sam = pysam.Samfile(local_sam_path, 'r')
    except:
        raise RuntimeError("[realignSamFile]Problem with SAM %s, exiting" % local_sam_path)

    # sort the cPecan results by batch, then discard the batch number. this is so they 'line up' with the sam
    alignment_hash      = getAlignedSegmentDictionary(sam)
    sorted_cPecan_fids  = sortResultsByBatch(cPecan_cigar_fileIds)
    temp_sam_filepath   = job.fileStore.getLocalTempFileName()
    output_sam_handle   = pysam.Samfile(temp_sam_filepath, 'wb', template=sam)
    failed_rebuilds     = 0
    successful_rebuilds = 0
    uid                 = uuid.uuid4().hex[:4]
    failed_alignments   = LocalFile(workdir=job.fileStore.getLocalTempDir(),
                                    filename="{lab}_failed_{uid}_realigned.bam"
                                             "".format(lab=config["sample_label"],
                                                       uid=uid))
    failed_sam_handle   = pysam.Samfile(failed_alignments.fullpathGetter(), 'wb', template=sam)

    for pA in cigar_iterator():
        ops = []
        aR  = alignment_hash[pA.contig2]
        if len(aR.cigar) > 0 and aR.cigar[0][0] == 5:
            # Add any hard clipped prefix
            ops.append(aR.cigar[0])
        if aR.query_alignment_start > 0:
            ops.append((4, aR.qstart))
        ops += map(lambda op : (op.type, op.length), pA.operationList)
        if aR.query_alignment_end < len(aR.query_sequence):
            ops.append((4, len(aR.query_sequence) - aR.query_alignment_end))
        if len(aR.cigar) > 1 and aR.cigar[-1][0] == 5:
            # Add any hard clipped suffix
            ops.append(aR.cigar[-1])

        # Checks the final operation list, correct for the read
        ops_len   = sum(map(lambda (type, length) : length if type in (0, 1, 4) else 0, ops))
        cig_len   = sum(map(lambda (type, length) : length if type in (0, 1, 4) else 0, aR.cigar))
        match_len = sum(map(lambda (type, length) : length if type in (0, 2) else 0, ops))
        ref_len   = aR.reference_end - aR.reference_start

        if ops_len != cig_len or match_len != ref_len:
            failed_rebuilds += 1
            aR.cigar = tuple(ops)
            failed_sam_handle.write(aR)
            continue

        successful_rebuilds += 1
        aR.cigar = tuple(ops)
        # Write out
        output_sam_handle.write(aR)

    job.fileStore.logToMaster("[rebuildSamJobFunction]{f} alignments failed {s} worked"
                              "".format(f=failed_rebuilds, s=successful_rebuilds))
    sam.close()
    output_sam_handle.close()
    failed_sam_handle.close()

    require(os.path.exists(temp_sam_filepath),
            "[rebuildSamJobFunction]out_sam_path does not exist at {}".format(temp_sam_filepath))
    require(os.path.exists(failed_alignments.fullpathGetter()),
            "[rebuildSamJobFunction]failed alignments does not exist at"
            " {}".format(failed_alignments.fullpathGetter()))
    if failed_rebuilds > 0:
        deliverOutput(job, failed_alignments, config["output_dir"])
    return job.fileStore.writeGlobalFile(temp_sam_filepath)
