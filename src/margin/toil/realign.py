"""Functions to realign a SAM file with the base-level HMM
"""
import sys
import os
import pysam
import uuid
import toil_lib.programs as tlp
from margin.toil.localFileManager import LocalFileManager
from margin.utils import samIterator, getExonerateCigarFormatString, getFastaDictionary, fastaWrite


DEBUG = True
DOCKER_DIR = "/data/"


class TempBatchFileHandler(object):
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
        self._handle  = open(self.path, 'wa')
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


def cPecanRealignJobFunction(job, global_config, job_config,
                             cPecan_image="quay.io/artrand/cpecanrealign"):
    # type: (toil.job.Job, dict<string, parameters>, dict<string, string>)
    """Runs Docker-ized cPecan Hmm
    """
    uid = uuid.uuid4().hex
    # make fasta with the config
    contig_fasta_fn = "ref.{}.fa".format(uid)
    contig_fasta    = job.fileStore.getLocalTempDir() + "/" + contig_fasta_fn

    with open(contig_fasta, 'w') as fH:
        fastaWrite(fH, job_config["contig_name"], job_config["contig_seq"])

    # get the files locally, we need:
    out_fid = job.fileStore.getLocalTempFile()
    fids_to_get = {  # map of FileStoreIDs to filenames
        job_config["exonerate_fid"]  : "exonerate_file.{}.tmp".format(uid),
        job_config["query_seqs_fid"] : "query_seqs.{}.tmp".format(uid),
        job_config["query_labs_fid"] : "query_labs.{}.tmp".format(uid),
        global_config["hmm_file"]    : "in_hmm.{}.tmp".format(uid),
    }
    localFiles = LocalFileManager(job, fids_to_get.keys(), fids_to_get)

    # run cPecan in a container
    exonerate_arg   = "--exonerate_file {}".format(DOCKER_DIR + localFiles.localFileName(job_config["exonerate_fid"]))
    reference_arg   = "--reference {} ".format(DOCKER_DIR + contig_fasta_fn)
    query_seqs_arg  = "--query_seqs {} ".format(DOCKER_DIR + localFiles.localFileName(job_config["query_seqs_fid"]))
    query_labs_arg  = "--query_labels {} ".format(DOCKER_DIR + localFiles.localFileName(job_config["query_labs_fid"]))
    hmm_arg         = "--hmm_file {} ".format(DOCKER_DIR + localFiles.localFileName(global_config["hmm_file"]))
    gap_gamma_arg   = "--gap_gamma {} ".format(global_config["gap_gamma"])
    match_gamma_arg = "--match_gamma {} ".format(global_config["match_gamma"])
    ## TODO left off here, need a good way to aggregate `tempOutputFiles`
    output_arg      = "--output_alignment_file {} ".format()

    #cPecan_parameters = []
    #tlp.docker_call(tool=cPecan_image, parameters=



def realignSamFile(job, config, chained_alignment_output):
    # type: (toil.job.Job, dict<string, string>, dict<string, string>)
    # get the sam file locally
    local_sam_path  = job.fileStore.readGlobalFile(chained_alignment_output["chained_alignment_FileStoreID"])
    reference_fasta = job.fileStore.readGlobalFile(config["reference_FileStoreID"])
    assert(os.path.exists(reference_fasta)), "[realignSamFile]ERROR was not able to download reference from FileStore"
    reference_map = getFastaDictionary(reference_fasta)

    try:
        sam = pysam.Samfile(local_sam_path, 'r')
    except:
        job.fileStore.logToMaster("[realignSamFile]Problem with SAM %s, exiting" % local_sam_path)
        # TODO maybe throw an exception or something? How does toil handle errors?
        return

    def send_alignment_batch():
        if exonerate_cigar_batch is not None:
            # import the files to the FileStore
            if DEBUG:
                job.fileStore.logToMaster("[realignSamFile::send_alignment_batch] adding child job")

            assert(exonerate_cigar_batch.exists()), "[realignSamFile]ERROR didn't find exonerate file"
            assert(query_seqs.exists()), "[realignSamFile]ERROR didn't find query seqs file"
            assert(query_labs.exists()), "[realignSamFile]ERROR didn't find query labels file"
            exonerate_cigar_batch.close()
            query_seqs.close()
            query_labs.close()
            exonerate_fid  = job.fileStore.writeGlobalFile(exonerate_cigar_batch.path_getter())
            query_seqs_fid = job.fileStore.writeGlobalFile(query_seqs.path_getter())
            query_labs_fid = job.fileStore.writeGlobalFile(query_labs.path_getter())
            job.fileStore.logToMaster("GOING TO ALIGN --> %s %s %s" % (exonerate_fid, query_seqs_fid, query_labs_fid))
            cPecan_config = {
                "exonerate_fid"  : exonerate_fid,
                "query_seqs_fid" : query_seqs_fid,
                "query_labs_fid" : query_labs_fid,
                "contig_seq"     : reference_map[contig_name],  # this might be very inefficient for large genomes..?
                "contig_name"    : contig_name,
            }
            job.addChildJobFn(cPecanRealignJobFunction, config, cPecan_config)
        else:  # mostly for initial conditions, do nothing
            pass

    total_seq_len         = sys.maxint  # send a batch when we have this many bases
    exonerate_cigar_batch = None        # send a batch of exonerate-formatted cigars
    query_seqs            = None        # file containing read sequences
    query_labs            = None        # file containing read labels (headers)
    contig_name           = None        # send a batch when we get to a new contig

    # this loop shards the sam and sends batches to be realigned
    for aligned_segment in samIterator(sam):
        if (total_seq_len > config["max_length_per_job"] or
           contig_name != sam.getrname(aligned_segment.reference_id)):
            # make a new batch of reads and cigars
            send_alignment_batch()  # send the previous batch to become a child job
            exonerate_cigar_batch = TempBatchFileHandler(job.fileStore.getLocalTempFile())
            query_seqs            = TempBatchFileHandler(job.fileStore.getLocalTempFile())
            query_labs            = TempBatchFileHandler(job.fileStore.getLocalTempFile())
            total_seq_len         = 0

        assert(exonerate_cigar_batch is not None), "[realignSamFile]ERROR exonerate_cigar_batch is NONE"
        assert(query_seqs is not None), "[realignSamFile]ERROR query_batch is NONE"
        assert(query_labs is not None), "[realignSamFile]ERROR query_labs is NONE"

        exonerate_cigar_batch.write(getExonerateCigarFormatString(aligned_segment, sam) + "\n")
        query_seqs.write(aligned_segment.query_sequence + "\n")
        query_labs.write(aligned_segment.query_name + "\n")
        # updates
        total_seq_len += len(aligned_segment.query_sequence)
        contig_name = sam.getrname(aligned_segment.reference_id)

    send_alignment_batch()

    job.fileStore.logToMaster("-->Follow on?")

    sam.close()
