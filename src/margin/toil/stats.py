"""Functions to generate summary statistics for an alignment
"""
from __future__ import print_function
import os
import numpy

import pandas as pd

from toil_lib import require

from margin.utils import ReadAlignmentStats

from localFileManager import LocalFile, deliverOutput
from alignment import splitLargeAlignment


def collectAlignmentStatsJobFunction(job, config, full_alignment_fid, output_label):
    smaller_alns, uid_to_read = splitLargeAlignment(job,
                                                    config["stats_alignment_batch_size"],
                                                    full_alignment_fid)
    disk = 2 * config["reference_FileStoreID"].size

    stat_shards = [job.addChildJobFn(marginStatsJobFunction,
                                     alignment.FileStoreID,
                                     uid_to_read,
                                     config["reference_FileStoreID"],
                                     config["sample_FileStoreID"],
                                     config["local_alignment"]).rv()
                   for alignment in smaller_alns]

    job.addFollowOnJobFn(deliverAlignmentStats, config, stat_shards, output_label, disk=disk)
    return


def marginStatsJobFunction(job, alignment_fid, uid_to_read, reference_fid, fastq_fid, local_alignment):
    # type (toil.job.Job, filestoreID, filestoreID, filestoreID, bool)
    def write_stats(ras):
        l = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (uid_to_read[ras.read_label],
                                                  str(ras.readIdentity()),
                                                  str(ras.alignmentIdentity()),
                                                  str(ras.readCoverage()),
                                                  str(ras.mismatchesPerAlignedBase()),
                                                  str(ras.deletionsPerReadBase()),
                                                  str(ras.insertionsPerReadBase()),
                                                  str(ras.readLength()))
        _handle.write(l)

    sam_file    = job.fileStore.readGlobalFile(alignment_fid)
    ref_file    = job.fileStore.readGlobalFile(reference_fid)
    reads_fastq = job.fileStore.readGlobalFile(fastq_fid)
    out_stats   = job.fileStore.getLocalTempFileName()
    aln_stats   = ReadAlignmentStats.getReadAlignmentStats(sam_file,
                                                           reads_fastq,
                                                           ref_file,
                                                           uid_to_read,
                                                           (not local_alignment))

    _handle = open(out_stats, "w")
    map(write_stats, aln_stats)
    _handle.close()
    require(os.path.exists(out_stats), "[marginStatsJobFunction2]Didn't make stats file")
    return job.fileStore.writeGlobalFile(out_stats)


def deliverAlignmentStats(job, config, stat_shards, output_label):
    def parse_stats(f):
        df = pd.read_table(f,
                           usecols=(0, 1, 2, 3, 4, 5, 6, 7),
                           names=["read_label", "readIdentity", "alignmentIdentity",
                                  "coverage", "mismatchesPerAlignedBase", "deletionsPerReadBase",
                                  "insertionsPerReadBase", "readLength"],
                           dtype={"read_label"               : numpy.str,
                                  "readIdentity"             : numpy.float64,
                                  "alignmentIdentity"        : numpy.float64,
                                  "coverage"                 : numpy.float64,
                                  "mismatchesPerAlignedBase" : numpy.float64,
                                  "deletionsPerReadBase"     : numpy.float64,
                                  "insertionsPerReadBase"    : numpy.float64,
                                  "readLength"               : numpy.int})
        return df

    stats    = pd.concat([parse_stats(f) for f in [job.fileStore.readGlobalFile(fid) for fid in stat_shards]])
    out_file = LocalFile(workdir=job.fileStore.getLocalTempDir(),
                         filename="{sample}_{out_label}_stats.txt"
                         "".format(sample=config["sample_label"], out_label=output_label))
    _handle  = open(out_file.fullpathGetter(), "w")
    header   = "#read\treadIdentity\talignmentIdentity\treadCoverage\tmismatchesPerAlignedBase\t"\
               "deletionsPerReadBase\tinsertionsPerReadBase\treadLength\n"
    line     = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
    _handle.write(header)
    for _, row in stats.iterrows():
        _handle.write(line % (row["read_label"],
                              row["readIdentity"],
                              row["alignmentIdentity"],
                              row["coverage"],
                              row["mismatchesPerAlignedBase"],
                              row["deletionsPerReadBase"],
                              row["insertionsPerReadBase"],
                              row["readLength"]))
    _handle.close()
    deliverOutput(job, out_file, config["output_dir"])
    return


def marginStatsJobFunctionOld(job, config, alignment_fid, output_label):
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
