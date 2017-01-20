import os
import cPickle

from itertools import islice

from toil_lib import require
from toil_lib.programs import docker_call

from margin.toil.realign import setupLocalFiles, DOCKER_DIR
from margin.toil.variantCall import VariantCalls
from margin.marginCallerLib import \
    loadHmmSubstitutionMatrix,\
    getNullSubstitutionMatrix,\
    calcBasePosteriorProbs,\
    vcfWrite,\
    vcfWriteHeader
from margin.utils import getFastaDictionary
from localFileManager import LocalFile, deliverOutput


def calculateAlignedPairsJobFunction(job, global_config, job_config, hmm, batch_number,
                                     cPecan_image="quay.io/artrand/cpecanrealign"):
    workdir, local_hmm, local_output, local_input_obj = setupLocalFiles(job, global_config, hmm)

    if global_config["debug"]:
        job.fileStore.logToMaster("[calculateAlignedPairsJobFunction]Getting aligned pairs "
                                  "for batch {num} for contig {contig} and for {nseqs} "
                                  "sequences".format(num=batch_number,
                                                     contig=job_config["contig_name"],
                                                     nseqs=len(job_config["query_labels"])))

    with open(local_input_obj.fullpathGetter(), "w") as fH:
        cPickle.dump(job_config, fH)

    # cPecan container flags:
    aP_arg     = "--alignedPairs"
    input_arg  = "--input={}".format(DOCKER_DIR + local_input_obj.filenameGetter())
    hmm_file   = "--hmm_file={}".format(DOCKER_DIR + local_hmm.filenameGetter())
    output_arg = "--output_posteriors={}".format(DOCKER_DIR + local_output.filenameGetter())
    margin_arg = "--no_margin"
    cPecan_params = [aP_arg, input_arg, hmm_file, output_arg]

    if global_config["no_margin"]:
        cPecan_params.append(margin_arg)

    try:
        docker_call(job=job, tool=cPecan_image, parameters=cPecan_params, work_dir=(workdir + "/"))
        if not os.path.exists(local_output.fullpathGetter()):
            return None
        result_fid = job.fileStore.writeGlobalFile(local_output.fullpathGetter())
        return result_fid
    except:
        return None


def callVariantsOnBatch(job, config, expectations_batch):
    BASES         = "ACGT"
    variant_calls = VariantCalls()
    contig_seqs   = getFastaDictionary(job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
    error_model   = loadHmmSubstitutionMatrix(job.fileStore.readGlobalFile(config["error_model_FileStoreID"]))
    evo_sub_mat   = getNullSubstitutionMatrix()

    for contig in expectations_batch:
        contig_expectations = expectations_batch[contig]
        for position in expectations_batch[contig]:
            ref_base   = contig_seqs[contig][position]  # the base at this position in the reference
            pos_probs  = contig_expectations[position]  # a list [P(A), P(C), P(G), P(T)]
            total_prob = sum(pos_probs)
            require(total_prob > 0, "[callVariantsWithAlignedPairsJobFunction]Total prob == 0")
            posterior_probs = calcBasePosteriorProbs(dict(zip(BASES, map(lambda x : float(x) / total_prob, pos_probs))),
                                                     ref_base, evo_sub_mat, error_model)
            for b in BASES:
                if b != ref_base and posterior_probs[b] >= config["variant_threshold"]:
                    variant_calls.Add(contig, position, b, posterior_probs[b])

    return variant_calls


def vcfWriteJobFunction(job, reference_fasta, contig_seqs, variant_calls):
    vcf_file = job.fileStore.getLocalTempFile()
    vcfWrite(reference_fasta, contig_seqs, variant_calls, vcf_file, write_header=False)
    return job.fileStore.writeGlobalFile(vcf_file)


def vcfWriteJobFunction2(job, config, variant_calls):
    contig_seqs = getFastaDictionary(job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
    vcf_file = job.fileStore.getLocalTempFile()
    vcfWrite(config["ref"],  # the URL to the reference we used
             contig_seqs,
             variant_calls,
             vcf_file,
             write_header=False)
    return job.fileStore.writeGlobalFile(vcf_file)


def combineVcfShardsJobFunction(job, config, vcf_fids, output_label):
    workdir    = job.fileStore.getLocalTempDir()
    output_vcf = LocalFile(workdir=workdir,
                           filename="{sample}_{out_label}.vcf".format(sample=config["sample_label"],
                                                                      out_label=output_label))
    vcfWriteHeader(output_vcf.fullpathGetter(), config["ref"])
    with open(output_vcf.fullpathGetter(), "a") as fH:
        for fid in vcf_fids:
            vcf_shard  = job.fileStore.readGlobalFile(fid)
            vcf_handle = open(vcf_shard, "r")
            for line in vcf_handle:
                fH.write(line)
            vcf_handle.close()
            job.fileStore.deleteLocalFile(fid)

    deliverOutput(job, output_vcf, config["output_dir"])


def writeAndDeliverVCF(job, config, nested_variant_calls, output_label):
    job.fileStore.logToMaster("[writeAndDeliverVCF]Starting final VCF output")
    contig_seqs = getFastaDictionary(job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
    vcf_shards = []
    for variant_call_batch in nested_variant_calls:
        vcf_shards.append(job.addChildJobFn(vcfWriteJobFunction, config["ref"], contig_seqs, variant_call_batch).rv())
    job.addFollowOnJobFn(combineVcfShardsJobFunction, config, vcf_shards, output_label)


# FIXME rename this function because it returns variant calls now
def marginalizePosteriorProbsJobFunction(job, config, alignment_shard, cPecan_alignedPairs_fids):
    # type(toil.job.Job, dict, FileStoreID, list<FileStoreId>) n.b. the alignment_shard is ignored
    """reads in the posteriors and marginalizes (reduces) over the columns of the alignment returns a
    python dict with the expectations at each position
    """
    posterior_fids = [x[0] for x in cPecan_alignedPairs_fids]
    posterior_fids = [x for x in posterior_fids if x is not None]
    job.fileStore.logToMaster("[marginalizePosteriorProbsJobFunction]Collating posteriors have {} files ..."
                              "".format(len(posterior_fids)))

    expectations_at_each_position = {}  # stores posterior probs
    for aP_fid in posterior_fids:
        posteriors_file = job.fileStore.readGlobalFile(aP_fid)
        fH = open(posteriors_file, "r")
        expectations = cPickle.load(fH)  # dict (key: contig string, value: expectations dict<int, list>)
        for reference in expectations:   # adding a loop here, but there will almost always be just 1 reference
            if reference not in expectations_at_each_position:
                expectations_at_each_position[reference] = {}
            for position in expectations[reference]:  # position is int, expectations[reference] is dict<int, list>
                if position >= alignment_shard.start and position < alignment_shard.end:
                    if position not in expectations_at_each_position[reference]:
                        expectations_at_each_position[reference][position] = [0.0, 0.0, 0.0, 0.0]
                    expectations_at_each_position[reference][position] = \
                        [x + y
                            for x, y, in
                            zip(expectations_at_each_position[reference][position], expectations[reference][position])]
        job.fileStore.deleteGlobalFile(posteriors_file)

    job.fileStore.logToMaster("[marginalizePosteriorProbsJobFunction]... done")

    variant_calls = job.addChildJobFn(callVariantsOnBatch, config, expectations_at_each_position).rv()
    #write_vcf_job = job.addFollowOnJobFn(vcfWriteJobFunction2, config, variant_calls)
    #batch_vcf     = write_vcf_job.rv()
    return variant_calls  # a VariantCalls object


def combinePosteriorProbsJobFunction(job, expectations):
    BASES = "ACGT"
    if len(expectations) == 1:
        return expectations[0]
    x = expectations[0]
    y = expectations[1]
    for k in y:
        if k not in x.keys():
            x[k] = dict(zip(BASES, [0.0] * len(BASES)))
        for b in BASES:
            x[k][b] += y[k][b]
    return x


def parallelReducePosteriorProbsJobFunction(job, expectation_batches):
    combined_expectations = []
    PAIR_STEP = 2  # break into chunks of two 
    for i in xrange(0, len(expectation_batches), PAIR_STEP):
        combined_expectations.append(job.addChildJobFn(combinePosteriorProbsJobFunction,
                                                       expectation_batches[i : i + 2]).rv())
    if len(combined_expectations) > 1:
        return job.addFollowOnJobFn(parallelReducePosteriorProbsJobFunction, combined_expectations).rv()
    else:
        return combined_expectations[0]


def callVariantsWithAlignedPairsJobFunction1(job, config, expectation_batches, output_label):
    job.fileStore.logToMaster("[callVariantsWithAlignedPairsJobFunction1]Reducing {} expectation batches..."
                              "".format(len(expectation_batches)))
    # combine the expectations. this is crappy code
    posterior_probs_map = {}
    for batch in expectation_batches:  # keys are contigs, values are dicts<position, [probs]>
        for reference in batch:
            if reference not in posterior_probs_map:
                posterior_probs_map[reference] = {}
            for position in batch[reference]:  # we are now iterating over the positions with posterior probs
                if position not in posterior_probs_map[reference]:
                    posterior_probs_map[reference][position] = [0.0, 0.0, 0.0, 0.0]
                posterior_probs_map[reference][position] = \
                    [x + y for x, y in zip(posterior_probs_map[reference][position], batch[reference][position])]
    job.fileStore.logToMaster("[callVariantsWithAlignedPairsJobFunction1]...Done")
    job.addFollowOnJobFn(callVariantsWithAlignedPairsJobFunction2,
                         config,
                         posterior_probs_map,
                         output_label)


def callVariantsWithAlignedPairsJobFunction2(job, config, expectations_at_each_position, output_label):
    def chunk_expectations():
        it = iter(expectations_at_each_position)
        for i in xrange(0, len(expectations_at_each_position),
                        config["max_variant_call_positions_per_job"]):
            yield {k: expectations_at_each_position[k]
                   for k in islice(it, config["max_variant_call_positions_per_job"])}

    variant_calls = []
    batch_number  = 0
    # XXX TODO need an outter loop here 
    # parallel map of calling variants 
    for chunk in chunk_expectations():
        variant_calls.append(job.addChildJobFn(callVariantsOnBatch, config, chunk).rv())
        batch_number += 1

    job.fileStore.logToMaster("[callVariantsWithAlignedPairsJobFunction]Made {} batches for variant calling"
                              "".format(batch_number))

    job.addFollowOnJobFn(writeAndDeliverVCF, config, variant_calls, output_label)
