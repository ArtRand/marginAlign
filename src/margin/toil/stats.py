"""Functions to generate summary statistics for an alignment
"""
from __future__ import print_function
import numpy
from margin.utils import ReadAlignmentStats
from localFileManager import LocalFile, deliverOutput


def marginStatsJobFunction(job, config, alignment_fid, output_label):
    def report(values, statisticName):
        if not config["noStats"]:
            print("Average" + statisticName, numpy.average(values), file=fH)
            print("Median" + statisticName, numpy.median(values), file=fH)
            print("Min" + statisticName, min(values), file=fH)
            print("Max" + statisticName, max(values), file=fH)
        if config["printValuePerReadAlignment"]:
            print("Values" + statisticName, "\t".join(map(str, values)), file=fH)

    # get the alignment
    sam_file    = job.fileStore.readGlobalFile(alignment_fid)
    reads_fastq = job.fileStore.readGlobalFile(config["sample_FileStoreID"])
    ref_file    = job.fileStore.readGlobalFile(config["reference_FileStoreID"])
    workdir     = job.fileStore.getLocalTempDir()
    out_stats   = LocalFile(workdir=workdir,
                            filename="{sample}_{out_label}_stats.txt".format(sample=config["sample_label"],
                                                                             out_label=output_label))
    readAlignmentStats = ReadAlignmentStats.getReadAlignmentStats(sam_file,
                                                                  reads_fastq,
                                                                  ref_file,
                                                                  globalAlignment=(not config["local_alignment"]))

    with open(out_stats.fullpathGetter(), 'a+') as fH:
        if config["identity"]:
            report(map(lambda rAS : rAS.identity(), readAlignmentStats), "Identity")

        if config["readCoverage"]:
            report(map(lambda rAS : rAS.readCoverage(), readAlignmentStats), "ReadCoverage")

        if config["mismatchesPerAlignedBase"]:
            report(map(lambda rAS : rAS.mismatchesPerAlignedBase(), readAlignmentStats), "MismatchesPerAlignedBase")

        if config["deletionsPerReadBase"]:
            report(map(lambda rAS : rAS.deletionsPerReadBase(), readAlignmentStats), "DeletionsPerReadBase")

        if config["insertionsPerReadBase"]:
            report(map(lambda rAS : rAS.insertionsPerReadBase(), readAlignmentStats), "InsertionsPerReadBase")

        if config["readLength"]:
            report(map(lambda rAS : rAS.readLength(), readAlignmentStats), "ReadLength")

    deliverOutput(job, out_stats, config["output_dir"])
