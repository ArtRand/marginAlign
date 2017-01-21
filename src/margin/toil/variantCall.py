import datetime
import numpy as np
import pandas as pd

from itertools import chain

from margin.utils import getFastaDictionary
from localFileManager import LocalFile, deliverOutput

VARIANT_CALL_COLS = ["contig", "position", "alt", "posterior_prob"]


class VariantCalls(object):
    def __init__(self, list_of_calls):
        self.data = pd.DataFrame(list_of_calls)

    def Add(self, contig, position, alt_base, posterior_prob):
        self.data.loc[len(self.data)] = (contig, position, alt_base, posterior_prob)

    def Put(self, i, contig, position, alt_base, posterior_prob):
        self.data.loc[i] = (contig, position, alt_base, posterior_prob)

    def Clean(self):
        self.data = self.data.dropna()


class VariantCall(object):
    def __init__(self, contig, position, alt, posterior_prob):
        self.data        = pd.DataFrame(columns=VARIANT_CALL_COLS, index=[0])
        self.data.loc[0] = (contig, position, alt, posterior_prob)


def concatVariantCalls(parent_job, list_of_variant_calls):
    calls = [pd.read_pickle(parent_job.fileStore.readGlobalFile(fid)) for fid in list_of_variant_calls]
    return pd.concat(calls).sort_values(["contig", "position"])


def print_header(handle, reference_label):
    handle.write("##fileformat=VCFv4.2\n")
    handle.write("##fileDate=" + str(datetime.datetime.now().date()).replace("-", "") + "\n")
    handle.write("##source=marginCaller\n")
    handle.write("##reference=" + reference_label + "\n")
    handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")


def makeVcfFromVariantCallsJobFunction2(job, config, all_variant_calls, output_label):
    def parse_table(variant_calls):
        df = pd.read_table(variant_calls,
                           usecols=(0, 1, 2, 3, 4),
                           names=["contig", "position", "ref", "alt", "probability"],
                           dtype={"contig"      : np.str,
                                  "position"    : np.int,
                                  "ref"         : np.str,
                                  "probability" : np.float64})
        return df

    def collate_region_calls():
        # download the files
        job.fileStore.logToMaster("[makeVcfFromVariantCallsJobFunction2]Downloading %s files"
                                  % len(all_variant_calls))
        local_files = [job.fileStore.readGlobalFile(fid) for fid in all_variant_calls]
        return pd.concat([parse_table(vc) for vc in local_files])

    variant_calls = collate_region_calls()
    result_file   = LocalFile(workdir=job.fileStore.getLocalTempDir(),
                              filename="{sample}_{out_label}.vcf".format(sample=config["sample_label"],
                                                                         out_label=output_label))
    fH = open(result_file.fullpathGetter(), "w")
    print_header(fH, config["ref"])
    for contig, contig_df in variant_calls.groupby(["contig"]):
        for pos, pos_df in contig_df.groupby(["position"]):
            alts     = ",".join(pos_df["alt"])
            prob     = ",".join([str(x) for x in pos_df["probability"]])
            fH.write("%s\t%s\t.\t%s\t%s\t.\tPASS\t%s\n" % (contig,
                                                           pos + 1,
                                                           pos_df["ref"].iloc[0],
                                                           alts,
                                                           prob))
    fH.close()
    deliverOutput(job, result_file, config["output_dir"])


def makeVcfFromVariantCalls(job, config, all_variant_calls, output_label):

    def print_line(contig, position, ref_base, alt_base, prob):
        record = "{contig}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{prob}\n".format(contig=contig,
                                                                              pos=int(position),
                                                                              ref=ref_base,
                                                                              alt=alt_base,
                                                                              prob=prob)
        fH.write(record)
        return

    result_file = LocalFile(workdir=job.fileStore.getLocalTempDir(),
                            filename="{sample}_{out_label}.vcf".format(sample=config["sample_label"],
                                                                       out_label=output_label))

    with open(result_file.fullpathGetter(), "w") as fH:
        print_header(fH)
        contig_hash = getFastaDictionary(job.fileStore.readGlobalFile(config["reference_FileStoreID"]))
        results     = concatVariantCalls(job, all_variant_calls)
        for contig, contig_df in results.groupby(["contig"]):
            for pos, pos_df in contig_df.groupby(["position"]):
                ref_base = contig_hash[contig][int(pos)]
                alts     = ",".join(pos_df["alt"])
                prob     = ",".join([str(x) for x in pos_df["posterior_prob"]])
                print_line(contig=contig, position=(pos + 1), ref_base=ref_base, alt_base=alts, prob=prob)

    deliverOutput(job, result_file, config["output_dir"])
