import cPickle
import os

from itertools import islice, chain

#from toil_lib import flatten
from toil_lib import require
from toil_lib.programs import docker_call

from margin.toil.realign import setupLocalFiles, DOCKER_DIR, sortResultsByBatch
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
                                  "for batch {num}".format(num=batch_number))

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


def writeAndDeliverVCF(job, config, variantCalls, output_label):
    contig_seqs   = getFastaDictionary(job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
    workdir       = job.fileStore.getLocalTempDir()
    output_vcf    = LocalFile(workdir=workdir,
                              filename="{sample}_{out_label}.vcf".format(sample=config["sample_label"],
                                                                         out_label=output_label))

    job.fileStore.logToMaster("[writeAndDeliverVCF]Got %s chunks" % len(variantCalls))
    job.fileStore.logToMaster("[writeAndDeliverVCF]One looks like %s" % variantCalls[0])
    variant_calls2 = chain.from_iterable(variantCalls)

    #job.fileStore.logToMaster("[writeAndDeliverVCF]Got %s variant calls" % len(variant_calls2))
    #job.fileStore.logToMaster("[writeAndDeliverVCF]One looks like %s" % variant_calls2[0])

    vcfWrite(config["ref"], contig_seqs, variant_calls2, output_vcf.fullpathGetter())
    require(os.path.exists(output_vcf.fullpathGetter()), "[callVariantsWithAlignedPairsJobFunction]Did not make temp VCF file")

    if config["debug"]:
        variant_calls3 = chain.from_iterable(variantCalls)
        vcf_calls = vcfRead(output_vcf.fullpathGetter())
        calls     = set(map(lambda x : (x[0], x[1] + 1, x[2]), list(variant_calls3)))
        require(vcf_calls == calls, "[callVariantsWithAlignedPairsJobFunction]vcf write error")

    deliverOutput(job, output_vcf, config["output_dir"])


def callVariantsWithAlignedPairsJobFunction(job, config, input_samfile_fid, output_label, cPecan_alignedPairs_fids):
    def chunk_expectations():
        it = iter(expectations_at_each_position)
        for i in xrange(0, len(expectations_at_each_position), vcf_workers):
            yield {k: expectations_at_each_position[k] for k in islice(it, vcf_workers)}

    BASES = "ACGT"
    sorted_alignedPair_fids = sortResultsByBatch(cPecan_alignedPairs_fids)

    if config["debug"]:
        job.fileStore.logToMaster("[callVariantsWithAlignedPairsJobFunction]Got {} sets of aligned pairs "
                                  "".format(len(cPecan_alignedPairs_fids)))

    expectations_at_each_position = {}  # stores posterior probs

    for aP_fid in sorted_alignedPair_fids:
        posteriors_file = job.fileStore.readGlobalFile(aP_fid)
        with open(posteriors_file, 'r') as fH:
            posteriors = cPickle.load(fH)
            for k in posteriors:
                if k not in expectations_at_each_position:
                    expectations_at_each_position[k] = dict(zip(BASES, [0.0] * len(BASES)))
                for b in BASES:
                    expectations_at_each_position[k][b] += posteriors[k][b]

    vcf_workers   = int(len(expectations_at_each_position) / 4)    
    variant_calls = []
    i = 0
    for chunk in chunk_expectations():
        job.fileStore.logToMaster("-->adding chunk %s" % i)
        variant_calls.append(job.addChildJobFn(callVariantsOnBatch, config, chunk).rv())
        i += 1
    #variant_calls = []
    #contig_seqs   = getFastaDictionary(job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
    #error_model   = loadHmmSubstitutionMatrix(job.fileStore.readGlobalFile(config["error_model_FileStoreID"]))
    #evo_sub_mat   = getNullSubstitutionMatrix()
    #workdir       = job.fileStore.getLocalTempDir()
    #output_vcf    = LocalFile(workdir=workdir,
    #                          filename="{sample}_{out_label}.vcf".format(sample=config["sample_label"],
    #                                                                     out_label=output_label))

    #chunks        = [expectations_at_each_position.iteritems()] * vcf_workers
    #generator     = (dict(ifilter(None, v)) for v in izip_longest(*chunks))
    #variant_calls = [job.addChildJobFn(callVariantsOnBatch, config, batch).rv() for batch in generator]

    job.addFollowOnJobFn(writeAndDeliverVCF, config, variant_calls, output_label)

    #for contig, position in expectations_at_each_position:
    #    ref_base     = contig_seqs[contig][position]
    #    expectations = expectations_at_each_position[(contig, position)]
    #    total_prob   = sum(expectations.values())
    #    require(total_prob > 0, "[callVariantsWithAlignedPairsJobFunction]Total prob == 0")
    #    posterior_probs = calcBasePosteriorProbs(dict(zip(BASES, map(lambda x : float(expectations[x]) / total_prob, BASES))),
    #                                             ref_base, evo_sub_mat, error_model)
    #    for b in BASES:
    #        if b != ref_base and posterior_probs[b] >= config["variant_threshold"]:
    #            variant_calls.append((contig, position, b, posterior_probs[b]))

    #vcfWrite(config["ref"], contig_seqs, variant_calls, output_vcf.fullpathGetter())
    #require(os.path.exists(output_vcf.fullpathGetter()), "[callVariantsWithAlignedPairsJobFunction]Did not make temp VCF file")

    #if config["debug"]:
    #    vcf_calls = vcfRead(output_vcf.fullpathGetter())
    #    calls     = set(map(lambda x : (x[0], x[1] + 1, x[2]), variant_calls))
    #    require(vcf_calls == calls, "[callVariantsWithAlignedPairsJobFunction]vcf write error")

    #deliverOutput(job, output_vcf, config["output_dir"])
