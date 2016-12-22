#!/usr/bin/env python
"""Wrapper script for cPecanRealign in Docker
"""
import os
import sys
import uuid
import cPickle
from argparse import ArgumentParser

print("running.. RUNNER!!")


def fasta_write(file_path, seq, seq_label):
    # type (string, string, string)
    chunk_size = 100
    with open(file_path, 'w') as fH:
        fH.write(">{header}\n".format(header=seq_label))
        for i in xrange(0, len(seq), chunk_size):
            fH.write("%s\n" % seq[i : i + chunk_size])


def main():
    parser = ArgumentParser()
    # input files
    parser.add_argument("--input", action="store", dest="input_params", required=True)
    parser.add_argument("--hmm_file", action="store", dest="hmm_file", required=True)
    # input parameters
    parser.add_argument("--gap_gamma", action="store", type=float, required=True)
    parser.add_argument("--match_gamma", action="store", type=float, required=True)
    # output and environment
    parser.add_argument("--output_alignment_file", action="store", dest="out", required=True)
    parser.add_argument("--work_dir", action="store", dest="work_dir", default="/data/", required=False)
    parser.add_argument("--test_container", action="store_true", dest="test_container", default=False)
    args = parser.parse_args()

    if args.test_container:
        os.system("cPecanRealign --help")
        sys.exit(0)
    with open(args.input_params, 'r') as fH:
        params = cPickle.load(fH)

    # intermediate files
    uid = uuid.uuid4().hex
    temp_read_fn      = args.work_dir + "read.{}.fa".format(uid)
    temp_reference_fn = args.work_dir + "reference.{}.fa".format(uid)
    exonerate_cigars  = params["exonerate_cigars"]
    read_sequences    = params["query_sequences"]
    read_labels       = params["query_labels"]

    # write the reference for cPecan 
    fasta_write(file_path=temp_reference_fn, seq=params["contig_seq"], seq_label=params["contig_name"])

    # this loop aligns each read, to the reference with the exonerate pairwise alignment as the
    # guide alignment
    for cig, seq, lab in zip(exonerate_cigars, read_sequences, read_labels):
        fasta_write(file_path=temp_read_fn, seq=seq, seq_label=lab)
        cmd = "echo \"{cig}\" | cPecanRealign {ref} {query} --diagonalExpansion=10 "\
              "--splitMatrixBiggerThanThis=3000 --loadHmm {hmm} --gapGamma={gapG} "\
              "--matchGamma={matchG} >> {out}"
        os.system(cmd.format(cig=cig, ref=temp_reference_fn, query=temp_read_fn, hmm=args.hmm_file,
                             gapG=args.gap_gamma, matchG=args.match_gamma, out=args.out))

if __name__ == "__main__":
    main()
