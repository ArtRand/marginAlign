"""Functions to realign a SAM file with the base-level HMM
"""
import sys
import os
import pysam
from margin.utils import samIterator, getExonerateCigarFormatString
from sonLib.bioio import fastaWrite


class TempBatchFile(object):
    """A struct containing the path and handle for a file, used to easily access the a file and
    it's contents
    """
    def __init__(self, path):
        self.path     = path
        self.is_open  = False
        self._handle  = None

    def handle_getter(self):
        if not self.is_open:
            self.__open()
        return self._handle

    def __open(self):
        self._handle  = open(self.path, 'w')
        self.is_open = True

    def path_getter(self):
        return self.path

    def close(self):
        if self.is_open:
            self._handle.close()

    def exists(self):
        return os.path.exists(self.path)

    def write(self, str):
        if not self.is_open:
            self.__open()
        self._handle.write(str)


def cPecanRealign_job(job, exonerate_cigar_fid, query_fid):
    job.fileStore.logToMaster("[cPecanRealign_job]Got exonerate %s and query %s" % (exonerate_cigar_fid, query_fid))
    # get the files locally


def realignSamFile(job, config, chained_alignment_output):
    # type: (toil.job.Job, dict<string, string>, dict<string, string>)
    job.fileStore.logToMaster("[realignSamFile]START")
    # get the sam file locally
    local_sam_path = job.fileStore.readGlobalFile(chained_alignment_output["chained_alignment_FileStoreID"])
    job.fileStore.logToMaster("[realignSamFile]Got local samfile %s" % local_sam_path)
    try:
        sam = pysam.Samfile(local_sam_path, 'r')
    except:
        job.fileStore.logToMaster("[realignSamFile]Problem with SAM %s, exiting" % local_sam_path)
        # TODO maybe throw an exception or something? How does toil handle errors?
        return

    def send_alignment_batch():
        if exonerate_cigar_batch is not None:
            # import the files to the FileStore
            assert(exonerate_cigar_batch.exists()), "[realignSamFile]ERROR didn't find exonerate file"
            assert(query_batch.exists()), "[realignSamFile]ERROR didn't find query batch file"
            exonerate_cigar_batch.close()
            query_batch.close()
            exonerate_fid = job.fileStore.writeGlobalFile(exonerate_cigar_batch.path_getter())
            query_fid     = job.fileStore.writeGlobalFile(query_batch.path_getter())
            job.fileStore.logToMaster("GOING TO ALIGN --> %s %s" % (exonerate_fid, query_fid))
            job.addChildJobFn(cPecanRealign_job, exonerate_fid, query_fid)
        else:  # mostly for initial conditions, do nothing
            pass

    total_seq_len         = sys.maxint  # send a batch when we have this many bases
    exonerate_cigar_batch = None        # send a batch of exonerate-formatted cigars
    query_batch           = None        # file containing FASTA of read sequences
    contig_name           = None        # send a batch when we get to a new contig

    # this loop shards the sam and sends batches to be realigned
    for aligned_segment in samIterator(sam):
        if (total_seq_len > config["max_length_per_job"] or
           contig_name != sam.getrname(aligned_segment.reference_id)):
            # make a new batch of reads and cigars
            send_alignment_batch()  # send the previous batch to become a child job
            exonerate_cigar_batch = TempBatchFile(job.fileStore.getLocalTempFile())
            query_batch           = TempBatchFile(job.fileStore.getLocalTempFile())
            total_seq_len         = 0

        assert(exonerate_cigar_batch is not None), "[realignSamFile]ERROR exonerate_cigar_batch is NONE"
        assert(query_batch is not None), "[realignSamFile]ERROR query_batch is NONE"

        exonerate_cigar_batch.write(getExonerateCigarFormatString(aligned_segment, sam) + "\n")
        fastaWrite(query_batch.handle_getter(), aligned_segment.query_name, aligned_segment.query_sequence)
        # updates
        total_seq_len += len(aligned_segment.query_sequence)
        contig_name = sam.getrname(aligned_segment.reference_id)

    send_alignment_batch()

    job.fileStore.logToMaster("-->Follow on?")

    sam.close()
