"""Functions for performing expectation maximization for a base-aligning
pair HMM
"""
import os
import sys
import pysam
import random
import uuid
import toil_lib.programs as tlp
from localFileManager import LocalFileManager, LocalFile
from cPecan.cPecanEm import Hmm
from margin.utils import getExonerateCigarFormatString, samIterator

DOCKER_DIR = "/data/"


def performBaumWelchOnSamJobFunction(job, config, chain_alignment_output):
    if config["hmm_file"] is not None or not config["random_start"]:  # normal EM
        job.addFollowOnJobFn(prepareBatchesJobFunction, config, chain_alignment_output)
    else:
        raise NotImplementedError
    return


def prepareBatchesJobFunction(job, config, chain_alignment_output):
    # type (toil.job.Job, dict<string, string>, dict<string, string>)
    """What does this function do?
    """

    def get_and_upload_model():
        # establishes the starting model, uploads it to the FileStore and returns the FileStoreID
        # returns the FileStoreID of the starting, working model
        if config["hmm_file"] is not None:
            # load the input model, normalize it, and make a copy that we use as the starting model
            # this way the input model is not changed 
            job.fileStore.logToMaster("[get_and_upload_model]Loading HMM from {}".format(config["hmm_file"]))
            local_hmm = job.fileStore.readGlobalFile(config["hmm_file"])
            assert(os.path.exists(local_hmm)), "[get_and_upload_model]ERROR couldn't find local HMM here {}"\
                                               "".format(local_hmm)

            hmm = Hmm.loadHmm(local_hmm)
            job.fileStore.logToMaster("[get_and_upload_model]Loaded model type {}".format(hmm.modelType))
            hmm.normalise()
        else:
            # no input model was provided, make a blank one
            assert(config["model_type"] is not None), "[get_and_upload_model]ERROR No model or model type provided"
            job.fileStore.logToMaster("[get_and_upload_model]Making model of type {}".format(config["model_type"]))
            hmm = Hmm(config["model_type"])
            if config["random_start"]:
                job.fileStore.logToMaster("[get_and_upload_model]Using random starting parameters")
                hmm.randomise()
            else:
                hmm.equalise()

        # TODO implement this
        if config["set_Jukes_Cantor_emissions"] is not None:
            hmm.setEmissionsToJukesCantor(float(config["set_Jukes_Cantor_emissions"]))

        starting_hmm = job.fileStore.getLocalTempFile()
        hmm.write(starting_hmm)
        return job.fileStore.writeGlobalFile(starting_hmm)

    def shard_alignments():
        # iterate over the aligned segments and make lists of tuples, each tuple should be
        # (exonerate cigar, query seq, alignment length) and the list should have a total 
        # cumulative length of `maxLengthPerjob` return the list of lists of tuples        
        cum_alignment_len = 0
        split_alignments   = []
        alignment_batch    = []  # list of aligned segments/aligned regions
        for aR in samIterator(sam):
            aln_length = aR.query_alignment_length
            alignment_batch.append(aR)
            cum_alignment_len += aln_length
            if cum_alignment_len >= config["max_length_per_job"]:  # make a batch and reset
                split_alignments.append((alignment_batch, cum_alignment_len))
                cum_alignment_len = 0
                alignment_batch   = []
        if alignment_batch != []:  # catch any remaining alignments
            split_alignments.append((alignment_batch, cum_alignment_len))

        return split_alignments  # list of tuples of lists [([aR...aR], aln_len)...]

    def sample_alignments(split_alignments):
        def exonerate_aligned_segments(aR_list, sam):
            exonerate_cigars = [getExonerateCigarFormatString(aR, sam) for aR in aR_list]
            return exonerate_cigars

        assert(isinstance(split_alignments, list)), "[sample_alignments]ERROR input must be a list"
        assert(config["max_sample_alignment_length"] > 0), "[sample_alignments]ERROR max alignment length to samples < 0"
        random.shuffle(split_alignments)
        # add batches to sampled_alignments until we reach the alotted alignment length
        sampled_alignments  = []
        total_sample_length = 0
        for aR_list, aln_len in split_alignments:
            if total_sample_length <= config["max_sample_alignment_length"]:
                sampled_alignments.append(exonerate_aligned_segments(aR_list, sam))
                total_sample_length += aln_len
            else:
                break
        assert(len(sampled_alignments) >= 1)
        job.fileStore.logToMaster("[sample_alignments]Sampled {total} alignment bases split into {batches} batches"
                                  "".format(total=total_sample_length, batches=len(sampled_alignments)))
        return sampled_alignments  # a list of lists of AlignmentShards

    def package_sampled_batches(sampled_alignments):
        def pack_up(batch):
            tmp = job.fileStore.getLocalTempFile()
            with open(tmp, 'w') as fH:
                for cigar in batch:
                    fH.write(cigar + "\n")
            return job.fileStore.writeGlobalFile(tmp)

        batch_fids = [pack_up(batch) for batch in sampled_alignments]
        assert(len(batch_fids) == len(sampled_alignments))
        return batch_fids  # a list of FileStoreIDs

    # handle the model
    working_model_fid = get_and_upload_model()
    # download the chained SAM, load it
    local_chained_sam = job.fileStore.readGlobalFile(chain_alignment_output["chained_alignment_FileStoreID"])
    assert(os.path.exists(local_chained_sam)), "[shard_alignments]ERROR didn't find local_chained_sam here "\
                                               "{}".format(local_chained_sam)
    sam = pysam.Samfile(local_chained_sam, 'r')
    aln_batch_fids = package_sampled_batches(sample_alignments(shard_alignments()))
    sam.close()
    job.addFollowOnJobFn(expectationMaximisationJobFunction, config, working_model_fid, aln_batch_fids)


def expectationMaximisationJobFunction(job, config, working_model_fid, aln_batch_fids, running_likelihood=None, iteration=0):
    if iteration < config["em_iterations"]:
        job.fileStore.logToMaster("[expectationMaximisationJobFunction]At iteration {}".format(iteration))
        expectations_fids = [
            job.addChildJobFn(getExpectationsJobFunction, aln_batch, config, working_model_fid)
            for aln_batch in aln_batch_fids]
        # do maximization next, advance iteration, add to running likelihood
    else:
        job.fileStore.logToMaster("DONE")


def getExpectationsJobFunction(job, batch_fid, config, working_model_fid,
                               cPecan_image="quay.io/artrand/cpecanrealign"):
    """This JobFunction runs the cPecan Docker container to collect expectations
    for a batch of alignments, it returns the FileStoreID for the expectations
    file"""
    # download the files we need to run the Docker
    fids_to_get = [
        batch_fid,                        # the batch of exonerate CIGARs
        config["reference_FileStoreID"],  # reference
        config["sample_FileStoreID"],     # reads
        working_model_fid,                # input hmm
    ]
    local_files = LocalFileManager(job=job, fileIds_to_get=fids_to_get)
    # make a temp file to use for the expectations
    uid = uuid.uuid4().hex
    expectations_file = LocalFile(workdir=local_files.workDir(), filename="expectations.{}.expectations".format(uid))

    # run the docker
    em_arg           = "--em"
    aln_arg          = "--aln_file={}".format(DOCKER_DIR + local_files.localFileName(batch_fid))
    reference_arg    = "--reference={}".format(DOCKER_DIR +
                                               local_files.localFileName(config["reference_FileStoreID"]))
    query_arg        = "--query={}".format(DOCKER_DIR + local_files.localFileName(config["sample_FileStoreID"]))
    hmm_arg          = "--hmm_file={}".format(DOCKER_DIR + local_files.localFileName(working_model_fid))
    expectations_arg = "--expectations={}".format(DOCKER_DIR + expectations_file.filenameGetter())
    cPecan_params    = [em_arg, aln_arg, reference_arg, query_arg, hmm_arg, expectations_arg]
    tlp.docker_call(tool=cPecan_image,
                    parameters=cPecan_params,
                    work_dir=local_files.workDir())

    # upload the file to the jobstore
    job.fileStore.logToMaster("[getExpectationsJobFunction]Finished DOCKER!")
    return job.fileStore.writeGlobalFile(expectations_file.fullpathGetter())
