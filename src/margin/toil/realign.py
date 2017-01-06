"""Functions to realign a SAM file with the base-level HMM
"""
import sys
import os
import pysam
import uuid
import cPickle
import toil_lib.programs as tlp
from toil_lib import require
from localFileManager import LocalFile
from margin.utils import samIterator, getExonerateCigarFormatString, getFastaDictionary
from sonLib.bioio import cigarRead

DEBUG = True
DOCKER_DIR = "/data/"


def setupLocalFiles(parent_job, global_config):
    # uid to be super sure we don't have any file collisions
    uid = uuid.uuid4().hex

    # need a local directory to pass to docker
    workdir = parent_job.fileStore.getLocalTempDir()
    # need the hmm file for cPecan and a local file for the output alignments
    local_hmm    = LocalFile(workdir=workdir, filename="hmm.{}.txt".format(uid))
    local_output = LocalFile(workdir=workdir, filename="cPecan_out.{}.txt".format(uid))

    # copy the hmm file from the FileStore to local workdir
    if global_config["EM"]:  # if we did EM, use the trained model
        assert(global_config["output_model"] is not None)  # double check
        hmm_model_fid = parent_job.fileStore.importFile(global_config["output_model"])
        assert(hmm_model_fid is not None), "[cPecanRealignJobFunction]ERROR importing trained model"
    else:
        hmm_model_fid = global_config["input_hmm_FileStoreID"]
        assert(hmm_model_fid is not None), "[cPecanRealignJobFunction]No input model and no EM"
    # read the HMM to the local file location
    parent_job.fileStore.readGlobalFile(fileStoreID=hmm_model_fid, userPath=local_hmm.fullpathGetter())

    local_input_obj = LocalFile(workdir=workdir, filename="cPecanInput.{}.pkl".format(uid))

    return workdir, local_hmm, local_output, hmm_model_fid, local_input_obj


def sortResultsByBatch(cPecan_result_fids):
    batch_sorted       = sorted(cPecan_result_fids, key=lambda tup: tup[1])
    sorted_cPecan_fids = [x[0] for x in batch_sorted]
    return sorted_cPecan_fids


def realignSamFileJobFunction(job, config, input_samfile_fid):
    job.addFollowOnJobFn(shardSamJobFunction, config, input_samfile_fid,
                         cPecanRealignJobFunction, rebuildSamJobFunction)


def cPecanRealignJobFunction(job, global_config, job_config, batch_number,
                             cPecan_image="quay.io/artrand/cpecanrealign"):
    # type: (toil.job.Job, dict<string, parameters>, dict<string, string>)
    """Runs Docker-ized cPecan HMM
    """
    # uid to be super sure we don't have any file collisions
    #uid = uuid.uuid4().hex

    # need a local directory to pass to docker
    #workdir = job.fileStore.getLocalTempDir()
    # need the hmm file for cPecan and a local file for the output alignments
    #local_hmm    = LocalFile(workdir=workdir, filename="hmm.{}.txt".format(uid))
    #local_output = LocalFile(workdir=workdir, filename="cPecan_out.{}.txt".format(uid))

    # copy the hmm file from the FileStore to local workdir
    #if global_config["EM"]:  # if we did EM, use the trained model
    #    assert(global_config["output_model"] is not None)  # double check
    #    hmm_model_fid = job.fileStore.importFile(global_config["output_model"])
    #    assert(hmm_model_fid is not None), "[cPecanRealignJobFunction]ERROR importing trained model"
    #else:
    #    hmm_model_fid = global_config["input_hmm_FileStoreID"]
    #    assert(hmm_model_fid is not None), "[cPecanRealignJobFunction]No input model and no EM"
    # read the HMM to the local file location
    #job.fileStore.readGlobalFile(fileStoreID=hmm_model_fid, userPath=local_hmm.fullpathGetter())
    workdir, local_hmm, local_output, hmm_model_fid, local_input_obj = setupLocalFiles(job, global_config)
    if DEBUG:
        job.fileStore.logToMaster("[cPecanRealignJobFunction]Batch {batch} using HMM from {fid} "
                                  "and EM is {em}".format(batch=batch_number,
                                                          fid=hmm_model_fid,
                                                          em=global_config["EM"]))


    # pickle the job_config, that contains the reference sequence, the query sequences, and 
    # the pairwise alignments in exonerate format
    #local_input_obj = LocalFile(workdir=workdir, filename="cPecanInput.{}.pkl".format(uid))
    with open(local_input_obj.fullpathGetter(), "w") as fH:
        cPickle.dump(job_config, fH)

    # run cPecan in a container
    input_arg         = "--input={}".format(DOCKER_DIR + local_input_obj.filenameGetter())
    hmm_arg           = "--hmm_file={}".format(DOCKER_DIR + local_hmm.filenameGetter())
    gap_gamma_arg     = "--gap_gamma={}".format(global_config["gap_gamma"])
    match_gamma_arg   = "--match_gamma={}".format(global_config["match_gamma"])
    output_arg        = "--output_alignment_file={}".format(DOCKER_DIR + local_output.filenameGetter())
    cPecan_parameters = [input_arg, hmm_arg, gap_gamma_arg, match_gamma_arg, output_arg]
    tlp.docker_call(tool=cPecan_image, parameters=cPecan_parameters, work_dir=(workdir + "/"))

    # import the result to the FileStore
    result_fid = job.fileStore.writeGlobalFile(local_output.fullpathGetter(), cleanup=False)
    return result_fid


def rebuildSamJobFunction(job, config, input_samfile_fid, cPecan_cigar_fileIds):
    if DEBUG:
        job.fileStore.logToMaster("[rebuildSamJobFunction]Rebuild chained SAM {chained} with alignments "
                                  "from {cPecan_fids}"
                                  "".format(chained=input_samfile_fid,
                                            cPecan_fids=cPecan_cigar_fileIds.__str__()))

    # iterates over the files, downloads them, and iterates over the alignments 
    def cigar_iterator():
        for fid in sorted_cPecan_fids:
            local_copy = job.fileStore.readGlobalFile(fid)
            assert(os.path.exists(local_copy)), "[cigar_iterator]ERROR couldn't find alignment {}".format(local_copy)
            for pA in cigarRead(open(local_copy, 'r')):
                yield pA
            job.fileStore.deleteLocalFile(fid)

    # download the chained sam
    local_sam_path = job.fileStore.readGlobalFile(input_samfile_fid)
    try:
        sam = pysam.Samfile(local_sam_path, 'r')
    except:
        job.fileStore.logToMaster("[realignSamFile]Problem with SAM %s, exiting" % local_sam_path)
        # TODO maybe throw an exception or something? How does toil handle errors?
        return
    # sort the cPecan results by batch, then discard the batch number. this is so they 'line up' with the sam
    sorted_cPecan_fids = sortResultsByBatch(cPecan_cigar_fileIds)
    #batch_sorted       = sorted(cPecan_cigar_fileIds, key=lambda tup: tup[1])
    #sorted_cPecan_fids = [x[0] for x in batch_sorted]

    output_sam_path = job.fileStore.getLocalTempFile()
    output_sam      = pysam.Samfile(output_sam_path, 'wh', template=sam)

    for aR, pA in zip(samIterator(sam), cigar_iterator()):
        ops = []
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
        assert sum(map(lambda (type, length) : length if type in (0, 1, 4) else 0, ops)) == \
            sum(map(lambda (type, length) : length if type in (0, 1, 4) else 0, aR.cigar))
        # Correct for the reference
        assert sum(map(lambda (type, length) : length if type in (0, 2) else 0, ops)) == \
            aR.reference_end - aR.reference_start

        aR.cigar = tuple(ops)
        # Write out
        output_sam.write(aR)

    if DEBUG:
        job.fileStore.logToMaster("[rebuildSamJobFunction]Finished rebuilding SAM, exporting to"
                                  "{}".format(config["output_sam_path"]))
    sam.close()
    output_sam.close()
    job.fileStore.exportFile(output_sam_path, config["output_sam_path"])


def shardSamJobFunction(job, config, input_samfile_fid, batch_job_function, followOn_job_function):
    # type: (toil.job.Job, dict<string, string>, string,
    #        JobFunctionWrappingJob, JobFunctionWrappingJob)
    # get the sam file locally
    local_sam_path  = job.fileStore.readGlobalFile(input_samfile_fid)
    reference_fasta = job.fileStore.readGlobalFile(config["reference_FileStoreID"])
    require(os.path.exists(reference_fasta),
            "[realignSamFile]ERROR was not able to download reference from FileStore")
    reference_map = getFastaDictionary(reference_fasta)

    try:
        sam = pysam.Samfile(local_sam_path, 'r')
    except:
        job.fileStore.logToMaster("[realignSamFile]Problem with SAM %s, exiting" % local_sam_path)
        # TODO maybe throw an exception or something? How does toil handle errors?
        return

    def send_alignment_batch(result_fids, batch_number):
        # type: (list<string>)
        # result_fids should be updated with the FileStoreIDs with the cPecan results
        if exonerate_cigar_batch is not None:
            if DEBUG:
                job.fileStore.logToMaster("[send_alignment_batch]Starting batch {batch} for contig {contig} "
                                          "".format(batch=batch_number, contig=contig_name))

            assert(len(exonerate_cigar_batch) == len(query_seqs)),\
                "[send_alignment_batch] len(exonerate_cigar_batch) != len(query_seqs)"
            assert(len(query_seqs) == len(query_labs)),\
                "[send_alignment_batch] len(query_seqs) != len(query_labs)"

            cPecan_config = {
                "exonerate_cigars" : exonerate_cigar_batch,
                "query_sequences"  : query_seqs,
                "query_labels"     : query_labs,
                "contig_seq"       : reference_map[contig_name],  # this might be very inefficient for large genomes..?
                "contig_name"      : contig_name,
            }
            result_id = job.addChildJobFn(batch_job_function, config, cPecan_config, batch_number).rv()
            result_fids.append((result_id, batch_number))
            return batch_number + 1
        else:  # mostly for initial conditions, do nothing
            return batch_number

    total_seq_len         = sys.maxint  # send a batch when we have this many bases
    exonerate_cigar_batch = None        # send a batch of exonerate-formatted cigars
    query_seqs            = None        # list containing read sequences
    query_labs            = None        # list containing read labels (headers)
    contig_name           = None        # send a batch when we get to a new contig
    cPecan_results        = []          # container with the FileStoreIDs of the re-alignment results
    batch_number          = 0           # ordering of the batches, so we can reassemble the new sam later

    # this loop shards the sam and sends batches to be realigned
    for aligned_segment in samIterator(sam):
        if (total_seq_len > config["max_length_per_job"] or
           contig_name != sam.getrname(aligned_segment.reference_id)):  # make a new batch of reads and cigars
            # send the previous batch to become a child job
            batch_number = send_alignment_batch(result_fids=cPecan_results, batch_number=batch_number)
            # start new batches
            exonerate_cigar_batch = []
            query_seqs            = []
            query_labs            = []
            total_seq_len         = 0

        assert(exonerate_cigar_batch is not None), "[realignSamFile]ERROR exonerate_cigar_batch is NONE"
        assert(query_seqs is not None), "[realignSamFile]ERROR query_batch is NONE"
        assert(query_labs is not None), "[realignSamFile]ERROR query_labs is NONE"

        exonerate_cigar_batch.append(getExonerateCigarFormatString(aligned_segment, sam) + "\n")
        query_seqs.append(aligned_segment.query_sequence + "\n")
        query_labs.append(aligned_segment.query_name + "\n")
        # updates
        total_seq_len += len(aligned_segment.query_sequence)
        contig_name = sam.getrname(aligned_segment.reference_id)

    send_alignment_batch(result_fids=cPecan_results, batch_number=batch_number)

    job.addFollowOnJobFn(followOn_job_function, config, input_samfile_fid, cPecan_results)

    sam.close()
