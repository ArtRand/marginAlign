#!/usr/bin/env python
"""Wrapper script for cPecanRealign in Docker
"""
import os
import uuid
import cPickle
from argparse import ArgumentParser

print("running.. RUNNER!! with pickle")


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
    parser.add_argument("--exonerate_file", action="store", dest="exonerate_file", required=True)
    parser.add_argument("--reference", action="store", dest="reference_fasta", required=True)
    parser.add_argument("--query_seqs", action="store", dest="query_seqs", required=True)
    parser.add_argument("--query_labels", action="store", dest="query_labels", required=True)
    parser.add_argument("--hmm_file", action="store", dest="hmm_file", required=True)
    # input parameters
    parser.add_argument("--gap_gamma", action="store", type=float, required=True)
    parser.add_argument("--match_gamma", action="store", type=float, required=True)
    # output and environment
    parser.add_argument("--output_alignment_file", action="store", dest="out", required=True)
    parser.add_argument("--work_dir", action="store", dest="work_dir", default="/data/", required=False)
    args = parser.parse_args()

    assert(os.path.exists(args.exonerate_file)), "[cPecan::runner]ERROR didn't find exonerate file "\
                                                 "looked {}".format(args.exonerate_file)
    assert(os.path.exists(args.reference_fasta)), "[cPecan::runner]ERROR didn't find reference FASTA "\
                                                  "looked {}".format(args.reference_fasta)
    assert(os.path.exists(args.query_file)), "[cPecan::runner]ERROR didn't find query reads file "\
                                             "looked {}".format(args.query_file)

    # intermediate files
    uid = uuid.uuidt().hex
    temp_read        = args.work_dir + "read{}.fa".format(uid)
    exonerate_cigars = open(args.exonerate_file, 'r')
    read_sequences   = open(args.query_seqs, 'r')
    read_labels      = open(args.query_labels, 'r')

    for cig, seq, lab in zip(exonerate_cigars, read_sequences, read_labels):
        fasta_write(file_path=temp_read, seq=seq, seq_label=lab)
        cmd = "echo \"{cig}\" | cPecanRealign {ref} {read} --diagonalExpansion=10 "\
              "--splitMatrixBiggerThanThis=3000 --loadHmm {hmm} --gapGamma={gapG} "\
              "--matchGamma={matchG} >> {out}"
        os.system(cmd.format(cig=cig, ref=args.reference_fasta, read=temp_read, hmm=args.hmm_file,
                             gapG=args.gap_gamma, matchG=args.match_gamma, out=args.out))

if __name__ == "__main__":
    main()
