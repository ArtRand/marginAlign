"""Functions for performing expectation maximization for a base-aligning
pair HMM
"""
import os
import sys
import pysam
import random
import uuid
from toil_lib.programs import docker_call
from toil_lib import require
from sonLib.bioio import fastaWrite
from localFileManager import LocalFileManager, LocalFile, deliverOutput
from hmm import Hmm, SYMBOL_NUMBER
from margin.utils import getExonerateCigarFormatString, samIterator

DOCKER_DIR = "/data/"


def performBaumWelchOnSamJobFunction(job, config, input_samfile_fid):
    disk   = int(2.5 * input_samfile_fid.size)
    memory = 6 * input_samfile_fid.size
    if config["input_hmm_FileStoreID"] is not None or not config["random_start"]:  # normal EM
        job.fileStore.logToMaster("[performBaumWelchOnSamJobFunction]Asking for disk {disk} and "
                                  "memory {mem} for batch prep".format(disk=disk, mem=memory))
        job.addFollowOnJobFn(prepareBatchesJobFunction, config, input_samfile_fid,
                             disk=disk, memory=memory)
    else:
        raise NotImplementedError
    return


def prepareBatchesJobFunction(job, config, input_samfile_fid):
    # type (toil.job.Job, dict<string, string>, dict<string, string>)
    """This JobFunction has three steps:
        1. Upload a model file to the FileStore to be the working model
           that we start with.
        2. Shard the alignments in the chained SAM into batches that have roughly
           'max_length_per_job' worth of alignment length to them.
        3. Sample from the alignment batches and
           prepare them in Exonerate/FASTA format for cPecan, upload the batches to the FileStore
           and send the fileIDs to the follow on function.
    """

    def get_and_upload_model():
        # establishes the starting model, uploads it to the FileStore and returns the FileStoreID
        # returns the FileStoreID of the starting, working model
        if config["input_hmm_FileStoreID"] is not None:
            # load the input model, normalize it, and make a copy that we use as the starting model
            # this way the input model is not changed 
            job.fileStore.logToMaster("[get_and_upload_model]Loading HMM from {}".format(config["input_hmm_FileStoreID"]))
            local_hmm = job.fileStore.readGlobalFile(config["input_hmm_FileStoreID"])
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
                job.fileStore.logToMaster("[get_and_upload_model]Using equal transition starting parameters")
                hmm.equalise()

        # TODO implement this
        if config["set_Jukes_Cantor_emissions"] is not None:
            hmm.setEmissionsToJukesCantor(float(config["set_Jukes_Cantor_emissions"]))

        starting_hmm = job.fileStore.getLocalTempFile()
        hmm.write(starting_hmm)
        return job.fileStore.writeGlobalFile(starting_hmm)

    def shard_alignments():
        # iterate over the aligned segments and make lists of tuples, each tuple should be
        # ([aligned_segments...], alignment_length)  
        cum_alignment_len = 0
        split_alignments  = []  # contains the batches
        alignment_batch   = []  # list of aligned segments/aligned regions
        for aR in samIterator(sam):
            aln_length = aR.query_alignment_length
            alignment_batch.append(aR)
            cum_alignment_len += aln_length
            if cum_alignment_len >= config["max_alignment_length_per_job"]:  # make a batch and reset
                split_alignments.append((alignment_batch, cum_alignment_len))
                cum_alignment_len = 0
                alignment_batch   = []
        if alignment_batch != []:  # catch any remaining alignments
            split_alignments.append((alignment_batch, cum_alignment_len))

        return split_alignments  # list of tuples of lists [([aR...aR], aln_len)...]

    def sample_alignments(split_alignments):
        assert(isinstance(split_alignments, list)), "[sample_alignments]ERROR input must be a list"
        assert(config["max_sample_alignment_length"] > 0), "[sample_alignments]ERROR max alignment length to samples < 0"
        random.shuffle(split_alignments)
        # add batches to sampled_alignments until we reach the alotted sample alignment length
        sampled_alignments  = []
        total_sample_length = 0
        for aR_list, aln_len in split_alignments:
            if total_sample_length <= config["max_sample_alignment_length"]:
                sampled_alignments.append(aR_list)
                total_sample_length += aln_len
            else:
                break
        assert(len(sampled_alignments) >= 1)
        job.fileStore.logToMaster("[sample_alignments]Sampled {total} alignment bases split into {batches} batches"
                                  "".format(total=total_sample_length, batches=len(sampled_alignments)))
        return sampled_alignments

    def package_sampled_batches(sampled_alignments):
        def pack_up(aR_list):
            cigar_file   = job.fileStore.getLocalTempFile()
            reads_file   = job.fileStore.getLocalTempFile()
            cigar_handle = open(cigar_file, 'w')
            reads_handle = open(reads_file, 'w')
            for aR in aR_list:
                cigar_handle.write(getExonerateCigarFormatString(aR, sam) + "\n")
                fastaWrite(reads_handle, aR.query_name, aR.seq)
            cigar_handle.close()
            reads_handle.close()
            cigar_fid = job.fileStore.writeGlobalFile(cigar_file)
            reads_fid = job.fileStore.writeGlobalFile(reads_file)
            # TODO test deleting local file here, to save on disk
            return (cigar_fid, reads_fid)

        batch_fids = [pack_up(batch) for batch in sampled_alignments]
        assert(len(batch_fids) == len(sampled_alignments))
        return batch_fids  # a list of tuples of FileStoreIDs (cigar_fid, fasta_fid)

    # handle the model
    working_model_fid = get_and_upload_model()
    # download the chained SAM, load it
    local_chained_sam = job.fileStore.readGlobalFile(input_samfile_fid)
    assert(os.path.exists(local_chained_sam)), "[shard_alignments]ERROR didn't find local_chained_sam here "\
                                               "{}".format(local_chained_sam)
    sam = pysam.Samfile(local_chained_sam, 'r')  # used throughout the nested functions
    # shard the alignments, sample alignments randomly, then prepare the formats and files for cPecan
    batch_fids = package_sampled_batches(sample_alignments(shard_alignments()))
    sam.close()
    job.addFollowOnJobFn(expectationMaximisationJobFunction, config, working_model_fid, batch_fids)


def expectationMaximisationJobFunction(job, config, working_model_fid, batch_fids,
                                       running_likelihood=None, iteration=0):
    if iteration < config["em_iterations"]:
        job.fileStore.logToMaster("[expectationMaximisationJobFunction]At iteration {}".format(iteration))
        expectations_fids = [
            job.addChildJobFn(getExpectationsJobFunction, batch, config, working_model_fid).rv()
            for batch in batch_fids]

        # do maximization next, advance iteration, add to running likelihood
        if running_likelihood is None:
            running_likelihood = []

        job.addFollowOnJobFn(maximizationJobFunction, config, expectations_fids, working_model_fid,
                             batch_fids, running_likelihood, iteration)
    else:
        job.fileStore.logToMaster("[expectationMaximisationJobFunction]Performed %s iterations" % iteration)
        # get local copy of the working_model
        if config["debug"]:
            job.fileStore.logToMaster("[expectationMaximisationJobFunction]Downloading trained model"
                                      " with FileStoreID {}".format(working_model_fid))
        trained_model_path = job.fileStore.readGlobalFile(working_model_fid, mutable=True)
        assert(os.path.exists(trained_model_path)), "[expectationMaximisationJobFunction]ERROR getting local "\
                                                    "copy of the trained model"
        # add the likelihoods to the bottom 
        with open(trained_model_path, 'a') as fH:
            fH.write("\t".join(map(str, running_likelihood)) + "\n")

        # upload the final, trained but unnormalized, model to the FileStore
        trained_model_fid = job.fileStore.writeGlobalFile(trained_model_path)
        if config["debug"]:
            job.fileStore.logToMaster("[expectationMaximisationJobFunction]Wrote final trained model to "
                                      "{}".format(trained_model_fid))
        job.fileStore.deleteGlobalFile(working_model_fid)
        config["unnormalized_model_FileStoreID"] = trained_model_fid
        job.addFollowOnJobFn(normalizeModelJobFunction, config)
        return


def getExpectationsJobFunction(job, batch_fid, config, working_model_fid,
                               cPecan_image="quay.io/artrand/cpecanrealign"):
    """This JobFunction runs the cPecan Docker container to collect expectations
    for a batch of alignments, it returns the FileStoreID for the expectations
    file"""
    # download the files we need to run the Docker
    assert(len(batch_fid) == 2), "[getExpectationsJobFunction]illegal batch_fid input"
    cigar_fid   = batch_fid[0]
    reads_fid   = batch_fid[1]
    fids_to_get = [
        cigar_fid,                        # CIGARS (in exonerate)
        reads_fid,                        # reads (in FASTA)
        config["reference_FileStoreID"],  # reference
        working_model_fid,                # input hmm
    ]
    local_files = LocalFileManager(job=job, fileIds_to_get=fids_to_get)
    # make a temp file to use for the expectations
    uid = uuid.uuid4().hex
    expectations_file = LocalFile(workdir=local_files.workDir(), filename="expectations.{}.expectations".format(uid))

    # run the docker
    em_arg           = "--em"
    aln_arg          = "--aln_file={}".format(DOCKER_DIR + local_files.localFileName(cigar_fid))
    reference_arg    = "--reference={}".format(DOCKER_DIR +
                                               local_files.localFileName(config["reference_FileStoreID"]))
    query_arg        = "--query={}".format(DOCKER_DIR + local_files.localFileName(reads_fid))
    hmm_arg          = "--hmm_file={}".format(DOCKER_DIR + local_files.localFileName(working_model_fid))
    expectations_arg = "--expectations={}".format(DOCKER_DIR + expectations_file.filenameGetter())
    cPecan_params    = [em_arg, aln_arg, reference_arg, query_arg, hmm_arg, expectations_arg]

    docker_call(job=job,
                tool=cPecan_image,
                parameters=cPecan_params,
                work_dir=local_files.workDir(),
                outfile=sys.stdout)

    # upload the file to the jobstore
    # XXX TODO XXX might need try/except guard here
    assert(os.path.exists(expectations_file.fullpathGetter())), "[getExpectationsJobFunction]Didn't find "\
                                                                "expectations file here"\
                                                                "{}".format(expectations_file.fullpathGetter())
    return job.fileStore.writeGlobalFile(expectations_file.fullpathGetter())


def maximizationJobFunction(job, config, expectations_fids, working_model_fid, aln_batch_fids,
                            running_likelihood, iteration):
    if len(expectations_fids) == 0:
        job.fileStore.logToMaster("[maximizationJobFunction]Didn't get any expectations FileStoreIDs")
        exit(1)

    if config["debug"]:
        job.fileStore.logToMaster("[maximizationJobFunction]Got %s expectations files" % len(expectations_fids))
        job.fileStore.logToMaster("[maximizationJobFunction]Based HMM on %s" % expectations_fids[0])

    # get the expetations files locally
    local_files = LocalFileManager(job, expectations_fids + [working_model_fid])
    hmm         = Hmm.loadHmm(local_files.localFilePath(expectations_fids[0]))
    for fid in expectations_fids[1:]:  # add them up and normalize
        hmm.addExpectationsFile(local_files.localFilePath(fid))
        job.fileStore.logToMaster("[maximizationJobFunction]Added %s" % fid)
    hmm.normalise()

    if config["debug"]:
        job.fileStore.logToMaster(
            "[maximizationJobFunction]On %i iteration got likelihood: %s for model-type: %s, model-file %s"
            % (iteration, hmm.likelihood, hmm.modelType, working_model_fid))
        job.fileStore.logToMaster(
            "[maximizationJobFunction]On %i iteration got transitions: %s for model-type: %s, model-file %s"
            % (iteration, " ".join(map(str, hmm.transitions)), hmm.modelType, working_model_fid))

    running_likelihood.append(hmm.likelihood)

    if config["train_emissions"]:
        hmm.tieEmissions()
        job.fileStore.logToMaster(
            "On %i iteration got emissions: %s for model-type: %s, model-file %s" %
            (iteration, " ".join(map(str, hmm.emissions)), hmm.modelType, working_model_fid))
    else:
        hmm.emissions = Hmm.loadHmm(local_files.localFilePath(working_model_fid)).emissions
        job.fileStore.logToMaster("On %i using the original emissions" % iteration)

    # TODO update band
    if config["update_band"]:
        raise NotImplementedError

    # write the new model
    new_model = job.fileStore.getLocalTempFileName()
    hmm.write(new_model)
    new_model_fid = job.fileStore.writeGlobalFile(new_model)
    job.fileStore.deleteGlobalFile(working_model_fid)
    job.addFollowOnJobFn(expectationMaximisationJobFunction,
                         config,
                         new_model_fid,
                         aln_batch_fids,
                         running_likelihood,
                         (iteration + 1))
    return


def normalizeModelJobFunction(job, config):
    def setHmmIndelEmissionsToBeFlat():
        """Set indel emissions to all be flat
        """
        for state in range(1, hmm.stateNumber):
            hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state + 1)] =\
                [1.0 / (SYMBOL_NUMBER**2)] * SYMBOL_NUMBER**2

    def normaliseHmmByReferenceGCContent(gcContent):
        """Normalise background emission frequencies to GC percent given
        """
        for state in range(hmm.stateNumber):
            # Don't normalise GC content of insert states (as they don't have any ref bases!)
            if state not in (2, 4):
                n = toMatrix(hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state + 1)])
                hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state + 1)] =\
                    fromMatrix(map(lambda i : map(lambda j : (n[i][j] / sum(n[i])) *
                               (gcContent / 2.0 if i in [1, 2] else (1.0 - gcContent) / 2.0),
                               range(SYMBOL_NUMBER)), range(SYMBOL_NUMBER)))  # Normalise

    assert("unnormalized_model_FileStoreID" in config.keys()),\
        "[normalizeModelJobFunction]unnormalized model FileStoreID not in config"

    if config["debug"]:
        job.fileStore.logToMaster("[normalizeModelJobFunction]normalizing model with FileStoreID %s"
                                  % config["unnormalized_model_FileStoreID"])

    toMatrix   = lambda e : map(lambda i : e[SYMBOL_NUMBER * i:SYMBOL_NUMBER * (i + 1)], xrange(SYMBOL_NUMBER))
    fromMatrix = lambda e : reduce(lambda x, y : list(x) + list(y), e)

    # get copy of unnormalized model
    unnormalized_model_path = job.fileStore.readGlobalFile(config["unnormalized_model_FileStoreID"])
    assert(os.path.exists(unnormalized_model_path)), "[normalizeModelJobFunction]ERROR getting model locally"
    hmm = Hmm.loadHmm(unnormalized_model_path)
    setHmmIndelEmissionsToBeFlat()
    normaliseHmmByReferenceGCContent(config["gc_content"])

    workdir = job.fileStore.getLocalTempDir()
    normalized_model = LocalFile(workdir=workdir, filename=Hmm.modelFilename(global_config=config))
    hmm.write(normalized_model.fullpathGetter())
    require(os.path.exists(normalized_model.fullpathGetter()),
            "[normalizeModelJobFunction]Didn't write VCF locally tried to write to {}"
            "".format(normalized_model.fullpathGetter()))

    deliverOutput(job, normalized_model, config["output_dir"])
    return
