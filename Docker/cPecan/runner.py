#!/usr/bin/env python
"""Wrapper script for cPecanRealign in Docker
"""
from __future__ import print_function
import os
import sys
import uuid
import cPickle
from argparse import ArgumentParser

print("running.. RUNNER!! with EM! and ALIGNEDPAIRS")


def fasta_write(file_path, seq, seq_label):
    # type (string, string, string)
    chunk_size = 100
    with open(file_path, 'w') as fH:
        fH.write(">{header}\n".format(header=seq_label))
        for i in xrange(0, len(seq), chunk_size):
            fH.write("%s\n" % seq[i : i + chunk_size])


def getParamsFromInputPickle(args):
    # check input
    assert(args.input_params is not None), "[runner.py::realign]No input arg (should be a pickle)"

    with open(args.input_params, 'r') as fH:
        params = cPickle.load(fH)

    return params


def getTempFiles(params, args):
    # intermediate files
    uid = uuid.uuid4().hex
    temp_read_fn      = args.work_dir + "read.{}.fa".format(uid)
    temp_reference_fn = args.work_dir + "reference.{}.fa".format(uid)
    exonerate_cigars  = params["exonerate_cigars"]
    read_sequences    = params["query_sequences"]
    read_labels       = params["query_labels"]
    return temp_read_fn, temp_reference_fn, exonerate_cigars, read_sequences, read_labels


def realign(args):
    # check input
    #assert(args.input_params is not None), "[runner.py::realign]No input arg (should be a pickle)"
    #with open(args.input_params, 'r') as fH:
    #    params = cPickle.load(fH)
    params = getParamsFromInputPickle(args)

    # intermediate files
    #uid = uuid.uuid4().hex
    #temp_read_fn      = args.work_dir + "read.{}.fa".format(uid)
    #temp_reference_fn = args.work_dir + "reference.{}.fa".format(uid)
    #exonerate_cigars  = params["exonerate_cigars"]
    #read_sequences    = params["query_sequences"]
    #read_labels       = params["query_labels"]
    temp_read_fn, temp_reference_fn, exonerate_cigars, read_sequences, read_labels = getTempFiles(params, args)

    # write the reference for cPecan 
    fasta_write(file_path=temp_reference_fn, seq=params["contig_seq"], seq_label=params["contig_name"])

    # this loop aligns each read, to the reference with the exonerate pairwise alignment as the
    # guide alignment
    for cig, seq, lab in zip(exonerate_cigars, read_sequences, read_labels):
        fasta_write(file_path=temp_read_fn, seq=seq, seq_label=lab)
        cmd = "echo \"{cig}\" | cPecanRealign {ref} {query} --diagonalExpansion=10 "\
              "--splitMatrixBiggerThanThis=3000 --loadHmm={hmm} --gapGamma={gapG} "\
              "--matchGamma={matchG} >> {out}"
        os.system(cmd.format(cig=cig, ref=temp_reference_fn, query=temp_read_fn, hmm=args.hmm_file,
                             gapG=args.gap_gamma, matchG=args.match_gamma, out=args.out))


def getAlignedPairs(args):
    BASES = "ACGT"

    params = getParamsFromInputPickle(args)

    temp_read_fn, temp_reference_fn, exonerate_cigars, read_sequences, read_labels = getTempFiles(params, args)

    fasta_write(file_path=temp_reference_fn, seq=params["contig_seq"], seq_label=params["contig_name"])

    # hash to store posterior probabilities in
    expectationsOfBasesAtEachPosition = {}

    uid = uuid.uuid4().hex
    temp_posterior_filepath = args.work_dir + "posteriorProbs.{}.txt".format(uid)
    for cig, seq, lab in zip(exonerate_cigars, read_sequences, read_labels):
        fasta_write(file_path=temp_read_fn, seq=seq, seq_label=lab)
        if args.no_margin:
            cmd = "echo \"{cig}\" | cPecanRealign {ref} {query} --diagonalExpansion=0 "\
                  "--splitMatrixBiggerThanThis=1 --rescoreOriginalAlignment "\
                  "--outputPosteriorProbs={probs}"
            os.system(cmd.format(cig=cig, ref=temp_reference_fn, query=temp_read_fn,
                                 probs=temp_posterior_filepath))
        else:
            cmd = "echo \"{cig}\" | cPecanRealign {ref} {query} --diagonalExpansion=10 "\
                  "--splitMatrixBiggerThanThis=100 --outputPosteriorProbs={probs} --loadHmm={hmm}"
            os.system(cmd.format(cig=cig, ref=temp_reference_fn, query=temp_read_fn,
                                 probs=temp_posterior_filepath, hmm=args.hmm_file))

        assert(os.path.exists(temp_posterior_filepath))

        # now collate the reference position expectations
        with open(temp_posterior_filepath, 'r') as fH:
            for refPosition, queryPosition, posteriorProb in map(lambda x : map(float, x.split()), fH):
                assert posteriorProb <= 1.01
                assert posteriorProb >= 0.0
                key = (params["contig_name"], int(refPosition))
                if key not in expectationsOfBasesAtEachPosition:
                    expectationsOfBasesAtEachPosition[key] = dict(zip(BASES, [0.0] * len(BASES)))
                queryBase = seq[int(queryPosition)].upper()
                if queryBase in BASES:  # Could be an N or other wildcard character, which we ignore
                    expectationsOfBasesAtEachPosition[key][queryBase] += 1.0 if args.no_margin else posteriorProb

    with open(args.out_posteriors, "w") as fH:
        cPickle.dump(expectationsOfBasesAtEachPosition, fH, cPickle.HIGHEST_PROTOCOL)


def em(args):
    assert(args.alignments is not None), "[runner.py::em]input alignments is None"
    assert(args.reference is not None), "[runner.py::em]input reference is None"
    assert(args.expectations is not None), "[runner.py::em]expectations outpath is None"
    assert(os.path.exists(args.hmm_file)), "[runner.py::em]Didn't find hmm file "\
                                           "file looked here {}".format(args.hmm_file)
    assert(os.path.exists(args.reference)), "[runner.py::em]Didn't find reference "\
                                            "file looked here {}".format(args.reference)
    assert(os.path.exists(args.query)), "[runner.py::em]Didn't find query "\
                                        "file looked here {}".format(args.query)
    assert(os.path.exists(args.alignments)), "[runner.py::em]Didn't find alignments "\
                                             "file looked here {}".format(args.alignments)

    cmd = "cat {alns} | cPecanRealign --logLevel DEBUG {reference} {reads} --loadHmm={hmm} "\
          "--outputExpectations={exps} --diagonalExpansion={diagEx} --splitMatrixBiggerThanThis={split}"\
          "".format(alns=args.alignments,
                    reads=args.query,
                    reference=args.reference,
                    hmm=args.hmm_file,
                    exps=args.expectations,
                    diagEx=args.diagonal_expansion,
                    split=args.split_matrix_threshold)

    print("[runner.py::em]Running {}".format(cmd), file=sys.stderr)

    os.system(cmd)

    return 0


def main():
    parser = ArgumentParser()
    # subprogram
    parser.add_argument("--em", action="store_true", dest="em", default=False)
    parser.add_argument("--alignedPairs", action="store_true", dest="alignedPairs", default=False)
    # input files
    parser.add_argument("--input", action="store", dest="input_params", required=False, default=None)
    parser.add_argument("--aln_file", action="store", dest="alignments", required=False, default=None)
    parser.add_argument("--query", action="store", dest="query", required=False, default=None)
    parser.add_argument("--reference", action="store", dest="reference", required=False, default=None)
    parser.add_argument("--expectations", action="store", dest="expectations", required=False, default=None)
    parser.add_argument("--hmm_file", action="store", dest="hmm_file", required=True)
    # input parameters
    parser.add_argument("--gap_gamma", action="store", type=float, required=False, default=0.5)
    parser.add_argument("--match_gamma", action="store", type=float, required=False, default=0.0)
    parser.add_argument("--diagonal_expansion", action="store", type=int, required=False, default=10)
    parser.add_argument("--split_matrix_threshold", action="store", type=int, required=False, default=300)
    parser.add_argument("--no_margin", action="store_true", dest="no_margin", required=False, default=False)
    # output and environment
    parser.add_argument("--output_alignment_file", action="store", dest="out", required=False,
                        default=None)
    parser.add_argument("--output_expectations", action="store", dest="expectations_out", required=False,
                        default=None)
    parser.add_argument("--output_posteriors", action="store", dest="out_posteriors", required=False,
                        default=None)
    parser.add_argument("--work_dir", action="store", dest="work_dir", default="/data/", required=False)
    args = parser.parse_args()

    if args.em:
        em(args)
    elif args.alignedPairs:
        getAlignedPairs(args)
    else:
        realign(args)


if __name__ == "__main__":
    main()
