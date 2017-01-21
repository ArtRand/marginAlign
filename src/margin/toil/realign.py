"""Functions to realign a SAM file with the base-level HMM
"""
import sys
import os
import pysam
import uuid
import cPickle
from toil_lib import require
from toil_lib.programs import docker_call
from alignment import splitLargeAlignment
from localFileManager import LocalFile, deliverOutput
from margin.toil.hmm import downloadHmm
from margin.utils import samIterator, getExonerateCigarFormatString, getFastaDictionary
from sonLib.bioio import cigarRead


DOCKER_DIR = "/data/"


def setupLocalFiles(parent_job, global_config, hmm):
    # uid to be super sure we don't have any file collisions
    uid             = uuid.uuid4().hex
    workdir         = parent_job.fileStore.getLocalTempDir()
    local_hmm       = LocalFile(workdir=workdir, filename="hmm.{}.txt".format(uid))
    local_output    = LocalFile(workdir=workdir, filename="cPecan_out.{}.txt".format(uid))
    local_input_obj = LocalFile(workdir=workdir, filename="cPecanInput.{}.pkl".format(uid))
    hmm.write(local_hmm.fullpathGetter())
    return workdir, local_hmm, local_output, local_input_obj


def sortResultsByBatch(cPecan_result_fids):
    batch_sorted       = sorted(cPecan_result_fids, key=lambda tup: tup[1])
    sorted_cPecan_fids = [x[0] for x in batch_sorted]
    return [x for x in sorted_cPecan_fids if x is not None]


def realignSamFileJobFunction(job, config, input_samfile_fid, output_label):
    smaller_alns        = splitLargeAlignment(job, config, input_samfile_fid)
    realigned_fids      = []
    hidden_markov_model = downloadHmm(job, config)

    for aln in smaller_alns:
        disk   = input_samfile_fid.size + config["reference_FileStoreID"].size
        memory = (6 * input_samfile_fid.size)
        realigned_sam = job.addChildJobFn(shardSamJobFunction, config, aln, hidden_markov_model,
                                          cPecanRealignJobFunction, rebuildSamJobFunction,
                                          disk=disk, memory=memory).rv()
        realigned_fids.append(realigned_sam)

    job.addFollowOnJobFn(combineRealignedSamfilesJobFunction, config, input_samfile_fid, realigned_fids,
                         output_label)


def combineRealignedSamfilesJobFunction(job, config, input_samfile_fid, realigned_fids, output_label):
    original_sam      = job.fileStore.readGlobalFile(input_samfile_fid)
    sam               = pysam.Samfile(original_sam, "r")
    filename          = "{sample}_{out_label}.bam".format(sample=config["sample_label"], out_label=output_label)
    output_sam        = LocalFile(workdir=job.fileStore.getLocalTempDir(), filename=filename)
    output_sam_handle = pysam.Samfile(output_sam.fullpathGetter(), "wb", template=sam)
    sam.close()

    for fid in realigned_fids:
        local_copy = job.fileStore.readGlobalFile(fid)
        samfile    = pysam.Samfile(local_copy, "rb")
        for line in samfile:
            output_sam_handle.write(line)
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
    with open(local_input_obj.fullpathGetter(), "w") as fH:
        cPickle.dump(job_config, fH)

    # run cPecan in a container
    input_arg         = "--input={}".format(DOCKER_DIR + local_input_obj.filenameGetter())
    hmm_arg           = "--hmm_file={}".format(DOCKER_DIR + local_hmm.filenameGetter())
    gap_gamma_arg     = "--gap_gamma={}".format(global_config["gap_gamma"])
    match_gamma_arg   = "--match_gamma={}".format(global_config["match_gamma"])
    output_arg        = "--output_alignment_file={}".format(DOCKER_DIR + local_output.filenameGetter())
    cPecan_parameters = [input_arg, hmm_arg, gap_gamma_arg, match_gamma_arg, output_arg]
    try:
        docker_call(job=job, tool=cPecan_image, parameters=cPecan_parameters, work_dir=(workdir + "/"))
        if not os.path.exists(local_output.fullpathGetter()):
            return None
        result_fid = job.fileStore.writeGlobalFile(local_output.fullpathGetter())
        return result_fid
    except:
        return None


def rebuildSamJobFunction(job, config, alignment_shard, cPecan_cigar_fileIds):
    alignment_fid = alignment_shard.FileStoreID
    if config["debug"]:
        job.fileStore.logToMaster("[rebuildSamJobFunction]Rebuild chained SAM {chained} with alignments "
                                  "from {cPecan_fids}"
                                  "".format(chained=alignment_fid,
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
    local_sam_path = job.fileStore.readGlobalFile(alignment_fid)
    try:
        sam = pysam.Samfile(local_sam_path, 'r')
    except:
        job.fileStore.logToMaster("[realignSamFile]Problem with SAM %s, exiting" % local_sam_path)
        # TODO maybe throw an exception or something? How does toil handle errors?
        return
    # sort the cPecan results by batch, then discard the batch number. this is so they 'line up' with the sam
    sorted_cPecan_fids = sortResultsByBatch(cPecan_cigar_fileIds)
    temp_sam_filepath  = job.fileStore.getLocalTempFileName()
    output_sam_handle  = pysam.Samfile(temp_sam_filepath, 'wb', template=sam)

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
        output_sam_handle.write(aR)

    sam.close()
    output_sam_handle.close()

    require(os.path.exists(temp_sam_filepath),
            "[rebuildSamJobFunction]out_sam_path does not exist at {}".format(temp_sam_filepath))
    return job.fileStore.writeGlobalFile(temp_sam_filepath)


def shardSamJobFunction(job, config, alignment_shard, hmm, batch_job_function, followOn_job_function):
    # get the sam file locally
    alignment_fid   = alignment_shard.FileStoreID
    local_sam_path  = job.fileStore.readGlobalFile(alignment_fid)
    reference_fasta = job.fileStore.readGlobalFile(config["reference_FileStoreID"])
    require(os.path.exists(reference_fasta),
            "[shardSamJobFunction]ERROR was not able to download reference from FileStore")
    reference_map = getFastaDictionary(reference_fasta)
    try:
        sam = pysam.Samfile(local_sam_path, 'r')
    except:
        job.fileStore.logToMaster("[realignSamFile]Problem with SAM %s, exiting" % local_sam_path)
        # TODO maybe throw an exception or something? How does toil handle errors?
        return

    def send_alignment_batch(result_fids, batch_number):
        # type: (list<string>, int) -> int
        # result_fids should be updated with the FileStoreIDs with the cPecan results
        if exonerate_cigar_batch is not None:
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

            # disk requirement = cigars + query_seqs + contig_seq + result
            # mem requirement = alignment
            result_id = job.addChildJobFn(batch_job_function, config, cPecan_config, hmm,
                                          batch_number, disk=config["reference_FileStoreID"].size).rv()
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
    alns_in_batch         = 0           # number of alignments we have in a batch, not to overload one

    # this loop shards the sam and sends batches to be realigned
    for aligned_segment in samIterator(sam):
        if (total_seq_len > config["max_alignment_length_per_job"] or
           contig_name != sam.getrname(aligned_segment.reference_id) or
           alns_in_batch >= config["max_alignments_per_job"] or
           len(aligned_segment.query_sequence) >= config["cut_batch_at_alignment_this_big"]):
            # send the previous batch to become a child job
            batch_number = send_alignment_batch(result_fids=cPecan_results, batch_number=batch_number)
            # start new batches
            exonerate_cigar_batch = []
            query_seqs            = []
            query_labs            = []
            total_seq_len         = 0
            alns_in_batch         = 0

        assert(exonerate_cigar_batch is not None), "[realignSamFile]ERROR exonerate_cigar_batch is NONE"
        assert(query_seqs is not None), "[realignSamFile]ERROR query_batch is NONE"
        assert(query_labs is not None), "[realignSamFile]ERROR query_labs is NONE"

        # TODO consider moving this to the child job, although you would have to send the sam too then...
        # unless I refactor the exonerateCigar function, which is totally doable, you just need the
        # reference name
        exonerate_cigar_batch.append(getExonerateCigarFormatString(aligned_segment, sam) + "\n")
        query_seqs.append(aligned_segment.query_sequence + "\n")
        query_labs.append(aligned_segment.query_name + "\n")
        # updates
        total_seq_len += len(aligned_segment.query_sequence)
        alns_in_batch += 1
        contig_name = sam.getrname(aligned_segment.reference_id)

    send_alignment_batch(result_fids=cPecan_results, batch_number=batch_number)
    sam.close()
    job.fileStore.logToMaster("[shardSamJobFunction]Made {} batches".format(batch_number + 1))
    # disk requirement <= alignment + exonerate cigars
    # memory requirement <= alignment 
    disk         = (2 * config["reference_FileStoreID"].size)
    memory       = (6 * alignment_fid.size)
    return job.addFollowOnJobFn(followOn_job_function, config, alignment_shard, cPecan_results, disk=disk,
                                memory=memory).rv()
