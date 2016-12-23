"""Functions for performing expectation maximization for a base-aligning
pair HMM
"""
import os
import pysam
import random
import cPickle
from cPecan.cPecanEm import Hmm
from margin.utils import getExonerateCigarFormatString, samIterator


class AlignmentShard(object):
    def __init__(self, exonerate_str, query_seq, aln_length, contig_label):
        self.exonerate_cigar = exonerate_str
        self.query_seq       = query_seq
        self.aln_length      = aln_length
        self.mapped_contig   = contig_label

    def exonerateCigar(self):
        return self.exonerate_cigar

    def querySeq(self):
        return self.query_seq

    def alnLength(self):
        return self.aln_length

    def mappedContig(self):
        return self.mapped_contig


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
        # download the chained SAM, load it
        local_chained_sam = job.fileStore.readGlobalFile(chain_alignment_output["chained_alignment_FileStoreID"])
        assert(os.path.exists(local_chained_sam)), "[shard_alignments]ERROR didn't find local_chained_sam here"\
                                                   "{}".format(local_chained_sam)
        sam = pysam.Samfile(local_chained_sam, 'r')

        # iterate over the aligned segments and make lists of tuples, each tuple should be
        # (exonerate cigar, query seq, alignment length) and the list should have a total 
        # cumulative length of `maxLengthPerjob` return the list of lists of tuples        
        cum_alignment_len = 0
        split_alignments   = []
        alignment_batch    = []
        for aR in samIterator(sam):
            exonerate_str = getExonerateCigarFormatString(aR, sam)
            query_seq     = aR.query_sequence
            aln_length    = aR.query_alignment_length
            contig_label  = sam.getrname(aR.reference_id)
            alignment_batch.append(AlignmentShard(exonerate_str, query_seq, aln_length, contig_label))
            #alignment_batch.append((exonerate_str, query_seq, aln_length, contig_label))
            cum_alignment_len += aln_length
            if cum_alignment_len >= config["max_length_per_job"]:
                split_alignments.append(alignment_batch)
                cum_alignment_len = 0
                alignment_batch = []
        if alignment_batch != []:
            split_alignments.append(alignment_batch)

        return split_alignments

    def sample_alignments(split_alignments):
        assert(isinstance(split_alignments, list)), "[sample_alignments]ERROR input must be a list"
        assert(config["max_sample_alignment_length"] > 0), "[sample_alignments]ERROR max alignment length to samples < 0"
        random.shuffle(split_alignments)
        # add batches to sampled_alignments until we reach the alotted alignment length
        sampled_alignments  = []
        total_sample_length = 0
        for aln_batch in split_alignments:
            if total_sample_length <= config["max_sample_alignment_length"]:
                sampled_alignments.append(aln_batch)
                total_sample_length += sum([shard.alnLength() for shard in aln_batch])
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
                for shard in batch:
                    fH.write(shard.exonerateCigar() + "\n")
            return job.fileStore.writeGlobalFile(tmp)

        batch_fids = [pack_up(batch) for batch in sampled_alignments]
        assert(len(batch_fids) == len(sampled_alignments))
        return batch_fids  # a list of FileStoreIDs

    working_model_fid = get_and_upload_model()
    aln_batch_fids    = package_sampled_batches(sample_alignments(shard_alignments()))

    job.addFollowOnJobFn(getExpectationsJobFunction, config, working_model_fid, aln_batch_fids)


def getExpectationsJobFunction(job, config, working_model_fid, aln_batch_fids, iteration=0):
    job.fileStore.logToMaster("GETTING EXPECTATIONS!!")

#def calcExpectations(job, )

