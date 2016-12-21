"""Functions to realign a SAM file with the base-level HMM
"""
import sys
import os
import pysam
import uuid
import cPickle
import toil_lib.programs as tlp
from margin.utils import samIterator, getExonerateCigarFormatString, getFastaDictionary


DEBUG = True
DOCKER_DIR = "/data/"


class LocalFile(object):
    """A struct containing the path and handle for a file, used to easily access the a file and
    it's contents
    """
    def __init__(self, workdir, filename):
        self.path     = workdir + "/" + filename
        self.filename = filename
        self.workdir  = workdir

    def filenameGetter(self):
        return self.filename

    def fullpathGetter(self):
        return self.path

    def workdirGetter(self):
        return self.workdir


def cPecanRealignJobFunction(job, global_config, job_config,
                             cPecan_image="quay.io/artrand/cpecanrealign"):
    # type: (toil.job.Job, dict<string, parameters>, dict<string, string>)
    """Runs Docker-ized cPecan Hmm
    """
    uid = uuid.uuid4().hex

    # need a local directory to pass to docker
    workdir = job.fileStore.getLocalTempDir()
    # need the hmm file for cPecan 
    local_hmm    = LocalFile(workdir=workdir, filename="hmm.{}.txt".format(uid))
    local_output = LocalFile(workdir=workdir, filename="cPecan_out.{}.txt".format(uid))
    # copy the hmm file from the FileStore to local workdir
    job.fileStore.readGlobalFile(fileStoreID=global_config["hmm_file"], userPath=local_hmm.fullpathGetter())
    # need a local file to put the output of the alignments

    if DEBUG:
        job.fileStore.logToMaster("local hmm %s local output %s" % (local_hmm.fullpathGetter(), local_output.fullpathGetter()))

    # pickle the job_config, that contains the reference sequence, the query sequences, and the pairwise alignments
    local_input_obj = LocalFile(workdir=workdir, filename="cPecanInput.{}.pkl".format(uid))
    with open(local_input_obj.fullpathGetter(), "w") as fH:
        cPickle.dump(job_config, fH)

    input_arg       = "--input={}".format(DOCKER_DIR + local_input_obj.filenameGetter())
    hmm_arg         = "--hmm_file={}".format(DOCKER_DIR + local_hmm.filenameGetter())
    gap_gamma_arg   = "--gap_gamma={}".format(global_config["gap_gamma"])
    match_gamma_arg = "--match_gamma={}".format(global_config["match_gamma"])
    output_arg      = "--output_alignment_file={}".format(DOCKER_DIR + local_output.filenameGetter())

    cPecan_parameters = [input_arg, hmm_arg, gap_gamma_arg, match_gamma_arg, output_arg]
    tlp.docker_call(tool=cPecan_image, parameters=cPecan_parameters, work_dir=(workdir + "/"))
    job.fileStore.logToMaster("DUMMY")


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
            if DEBUG:
                job.fileStore.logToMaster("[realignSamFile::send_alignment_batch] adding child job")

            # import the files to the FileStore
            assert(len(exonerate_cigar_batch) == len(query_seqs)), "[send_alignment_batch]"\
                " len(exonerate_cigar_batch) != len(query_seqs)"
            assert(len(query_seqs) == len(query_labs)), "[send_alignment_batch]"\
                " len(query_seqs) != len(query_labs)"
            job.fileStore.logToMaster("GOING TO ALIGN --> %s " % contig_name)
            cPecan_config = {
                "exonerate_cigars" : exonerate_cigar_batch,
                "query_sequences"  : query_seqs,
                "query_labels"     : query_labs,
                "contig_seq"       : reference_map[contig_name],  # this might be very inefficient for large genomes..?
                "contig_name"      : contig_name,
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

    send_alignment_batch()

    job.fileStore.logToMaster("-->Follow on?")

    sam.close()
