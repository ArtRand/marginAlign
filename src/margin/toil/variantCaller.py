import cPickle
import os

from itertools import islice, chain

from toil_lib import require
from toil_lib.programs import docker_call

from margin.toil.realign import setupLocalFiles, DOCKER_DIR
from margin.marginCallerLib import \
    loadHmmSubstitutionMatrix, getNullSubstitutionMatrix, calcBasePosteriorProbs, vcfWrite, vcfRead
from margin.utils import getFastaDictionary
from localFileManager import LocalFile, deliverOutput


def calculateAlignedPairsJobFunction(job, global_config, job_config, batch_number,
                                     cPecan_image="quay.io/artrand/cpecanrealign"):
    # UID to avoid file collisions
    workdir, local_hmm, local_output, local_input_obj = setupLocalFiles(job, global_config)

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

    docker_call(tool=cPecan_image, parameters=cPecan_params, work_dir=(workdir + "/"))

    require(os.path.exists(local_output.fullpathGetter()), "[calculateAlignedPairsJobFunction] "
            "Cannot find output aligned pairs")

    result_fid = job.fileStore.writeGlobalFile(local_output.fullpathGetter(), cleanup=False)
    return result_fid


def callVariantsOnBatch(job, config, expectations_batch):
    job.fileStore.logToMaster("[callVariantsOnBatch]working on a batch of expectations")
    BASES         = "ACGT"
    variant_calls = []
    contig_seqs   = getFastaDictionary(job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
    error_model   = loadHmmSubstitutionMatrix(job.fileStore.readGlobalFile(config["error_model_FileStoreID"]))
    evo_sub_mat   = getNullSubstitutionMatrix()

    for contig, position in expectations_batch:
        ref_base     = contig_seqs[contig][position]
        expectations = expectations_batch[(contig, position)]
        total_prob   = sum(expectations.values())
        require(total_prob > 0, "[callVariantsWithAlignedPairsJobFunction]Total prob == 0")
        posterior_probs = calcBasePosteriorProbs(dict(zip(BASES, map(lambda x : float(expectations[x]) / total_prob, BASES))),
                                                 ref_base, evo_sub_mat, error_model)
        for b in BASES:
            if b != ref_base and posterior_probs[b] >= config["variant_threshold"]:
                variant_calls.append((contig, position, b, posterior_probs[b]))

    return variant_calls


def writeAndDeliverVCF(job, config, nested_variant_calls, output_label):
    job.fileStore.logToMaster("[writeAndDeliverVCF]Starting final VCF output")
    contig_seqs = getFastaDictionary(job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
    workdir     = job.fileStore.getLocalTempDir()
    output_vcf  = LocalFile(workdir=workdir,
                            filename="{sample}_{out_label}.vcf".format(sample=config["sample_label"],
                                                                       out_label=output_label))

    variant_calls = chain.from_iterable(nested_variant_calls)

    vcfWrite(config["ref"], contig_seqs, variant_calls, output_vcf.fullpathGetter())
    require(os.path.exists(output_vcf.fullpathGetter()),
            "[callVariantsWithAlignedPairsJobFunction]Did not make temp VCF file")

    if config["debug"]:
        variant_calls2 = chain.from_iterable(nested_variant_calls)
        vcf_calls = vcfRead(output_vcf.fullpathGetter())
        calls     = set(map(lambda x : (x[0], x[1] + 1, x[2]), list(variant_calls2)))
        require(vcf_calls == calls, "[callVariantsWithAlignedPairsJobFunction]vcf write error")

    deliverOutput(job, output_vcf, config["output_dir"])


def marginalizePosteriorProbsJobFunction(job, config, input_samfile_fid, cPecan_alignedPairs_fids):
    # type(toil.job.Job, dict, FileStoreID, list<FileStoreId>) n.b. the input_samfile_fid is ignored
    """reads in the posteriors and marginalizes (reduces) over the columns of the alignment returns a 
    python dict with the expectations at each position
    """
    BASES = "ACGT"
    posterior_fids = [x[0] for x in cPecan_alignedPairs_fids]
    job.fileStore.logToMaster("[marginalizePosteriorProbsJobFunction]Collating posteriors have {} files ..."
                              "".format(len(posterior_fids)))

    expectations_at_each_position = {}  # stores posterior probs
    # run through the alignedPairs and collect the expectations for each (contig_name, position) key 
    # in this dict, this could be sped up hugely, but I'm not sure if it needs to be right now
    for aP_fid in posterior_fids:
        posteriors_file = job.fileStore.readGlobalFile(aP_fid)
        with open(posteriors_file, 'r') as fH:
            posteriors = cPickle.load(fH)
            for k in posteriors:
                if k not in expectations_at_each_position:
                    expectations_at_each_position[k] = dict(zip(BASES, [0.0] * len(BASES)))
                for b in BASES:
                    expectations_at_each_position[k][b] += posteriors[k][b]
        job.fileStore.deleteGlobalFile(posteriors_file)

    job.fileStore.logToMaster("[marginalizePosteriorProbsJobFunction]... done")

    return expectations_at_each_position


def callVariantsWithAlignedPairsJobFunction(job, config, input_samfile_fid, output_label,
                                            expectation_batches):
    def chunk_expectations():
        it = iter(expectations_at_each_position)
        for i in xrange(0, len(expectations_at_each_position),
                        config["max_variant_call_positions_per_job"]):
            yield {k: expectations_at_each_position[k]
                   for k in islice(it, config["max_variant_call_positions_per_job"])}

    job.fileStore.logToMaster("[callVariantsWithAlignedPairsJobFunction]Reducing {} expectations..."
                              "".format(len(expectation_batches)))
    BASES = "ACGT"
    # combine the expectations. this is crappy code
    expectations_at_each_position = expectation_batches[0]  # get the first one to reduce on
    for batch in expectation_batches[1:]:
        for k in batch:  # k = (contig, position)
            if k not in expectations_at_each_position:
                expectations_at_each_position[k] = batch[k]
                continue
            for b in BASES:
                expectations_at_each_position[k][b] += batch[k][b]

    job.fileStore.logToMaster("[callVariantsWithAlignedPairsJobFunction]...Done")
    variant_calls = []
    batch_number  = 0

    # parallel map of calling variants 
    for chunk in chunk_expectations():
        variant_calls.append(job.addChildJobFn(callVariantsOnBatch, config, chunk).rv())
        batch_number += 1

    job.fileStore.logToMaster("[callVariantsWithAlignedPairsJobFunction]Made {} batches for variant calling"
                              "".format(batch_number))

    job.addFollowOnJobFn(writeAndDeliverVCF, config, variant_calls, output_label)
