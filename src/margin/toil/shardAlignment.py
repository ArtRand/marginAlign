import sys
import os
import pysam

from bd2k.util.humanize import human2bytes
from toil_lib import require

from margin.utils import getFastaDictionary, samIterator, getExonerateCigarFormatStringWithCheck


def shardSamJobFunction(job, config, alignment_shard, hmm, batch_job_function, followOn_job_function,
                        exonerateCigarStringFn=getExonerateCigarFormatStringWithCheck,
                        batch_disk=human2bytes("1G"), followOn_disk=human2bytes("3G"),
                        batch_mem=human2bytes("2G"), followOn_mem=human2bytes("2G")):
    # type (toil.job.Job, dict, AlignmentShard, Hmm, JobFunction, JobFuncton, bytes, bytes, bytes, bytes)
    """workhorse job function that spawns alignment jobs. given an alignment shard (contains start, end
    and fileStore_ID) this function spaws off `batch_job_functions`s that perform work. it collects
    the return values from the `batch_job_function`s and hands them off to the followOn_job_function
    returns: the return value from the followOn_function (.rv())
    """
    alignment_fid   = alignment_shard.FileStoreID
    local_sam_path  = job.fileStore.readGlobalFile(alignment_fid)
    reference_fasta = job.fileStore.readGlobalFile(config["reference_FileStoreID"])
    require(os.path.exists(reference_fasta),
            "[shardSamJobFunction]ERROR was not able to download reference from FileStore")
    reference_map = getFastaDictionary(reference_fasta)
    try:
        sam = pysam.Samfile(local_sam_path, 'r')
    except:
        raise RuntimeError("[shardSamJobFunction]Problem opening alignment %s" % local_sam_path)

    def send_alignment_batch(result_fids, batch_number):
        # type: (list<string>, int) -> int
        # result_fids should be updated with the FileStoreIDs with the cPecan results
        if exonerate_cigar_batch is not None:
            assert(len(exonerate_cigar_batch) == len(query_seqs)),\
                "[send_alignment_batch] len(exonerate_cigar_batch) != len(query_seqs)"
            assert(len(query_seqs) == len(query_labs)),\
                "[send_alignment_batch] len(query_seqs) != len(query_labs)"
            # TODO add try/except for contigs that aren't in the reference map
            cPecan_config = {
                "exonerate_cigars" : exonerate_cigar_batch,
                "query_sequences"  : query_seqs,
                "query_labels"     : query_labs,
                "contig_seq"       : reference_map[contig_name],  # this might be very inefficient for large genomes..?
                "contig_name"      : contig_name,
            }

            result_id = job.addChildJobFn(batch_job_function, config, cPecan_config, hmm,
                                          batch_number, disk=batch_disk, memory=batch_mem).rv()
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

        #exonerate_cigar_batch.append(getExonerateCigarFormatString(aligned_segment, sam) + "\n")
        exonerate_cigar, ok = exonerateCigarStringFn(aligned_segment, sam)
        if not ok:
            continue
        #exonerate_cigar_batch.append(exonerateCigarStringFn(aligned_segment, sam) + "\n")
        exonerate_cigar_batch.append(exonerate_cigar + "\n")
        query_seqs.append(aligned_segment.query_sequence + "\n")
        query_labs.append(aligned_segment.query_name + "\n")
        # updates
        total_seq_len += len(aligned_segment.query_sequence)
        alns_in_batch += 1
        contig_name = sam.getrname(aligned_segment.reference_id)

    send_alignment_batch(result_fids=cPecan_results, batch_number=batch_number)
    sam.close()
    job.fileStore.logToMaster("[shardSamJobFunction]Made {} batches".format(batch_number + 1))
    return job.addFollowOnJobFn(followOn_job_function, config, alignment_shard, cPecan_results,
                                disk=followOn_disk, memory=followOn_mem).rv()
