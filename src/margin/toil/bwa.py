#!/usr/bin/env python
"""Burrow-Wheeler Aligner for pairwise alignment between DNA sequences
Module for running a BWA alignment program in a docker container.
Cite: Li H. (2013) Aligning sequence reads, clone sequences and assembly
contigs with BWA-MEM. arXiv:1303.3997v2
"""
from __future__ import print_function
import os
import uuid
from itertools import izip
import toil_lib.programs as tlp
from localFileManager import LocalFileManager
from alignment import AlignmentStruct, AlignmentFormat

DEBUG = True
DOCKER_DIR = "/data/"


def bwa_index_file_suffixes():
    return [".amb", ".ann", ".bwt", ".pac", ".sa"]


def bwa_index_docker_call(job, bwa_fileId_map,
                          memory="10M", cores=1, disk="10M",
                          bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    # type: (toil.job.Job, string, string, 
    #        bool, 
    #        string, int, string, 
    #        string)
    """Builds all required indices (crufty files) that BWA needs and
    imports them into the fileStore
    """
    def _run_bwa_index():
        bwa_index_parameters = ["index", dkr_reference_path]
        if DEBUG:
            job.fileStore.logToMaster("[bwa_index_docker_call::_run_bwa_index]workDir: {}".format(localFiles.workDir()))

        tlp.docker_call(tool=bwa_docker_image,
                        parameters=bwa_index_parameters,
                        work_dir=localFiles.workDir())

        # bwa docker call creates index files in the local working directory
        # import files to fileStore
        bwa_index_urls = [localFiles.localFilePath(bwa_fileId_map["reference_fasta"]) +
                          suffix for suffix in bwa_index_file_suffixes()]
        new_ids        = [job.fileStore.writeGlobalFile(x) for x in bwa_index_urls]

        # make a map of suffix files to their file Id
        index_fileId_map = dict([(suf, fid) for suf, fid in izip(bwa_index_file_suffixes(), new_ids)])

        return index_fileId_map

    # get a local copy of the reference and reads files for docker
    localFiles = LocalFileManager(job, [bwa_fileId_map["reference_fasta"]])

    # arguments for the bwa indexing and alignment
    dkr_reference_path = DOCKER_DIR + localFiles.localFileName(bwa_fileId_map["reference_fasta"])
    return _run_bwa_index()


def bwa_docker_align(job, bwa_input_map, bwa_index_map, bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    # move the read and reference files to local working directory
    uid1, uid2, uid3 = uuid.uuid4().hex, uuid.uuid4().hex, uuid.uuid4().hex
    user_paths = {  # map of fileIds to desired file names
        bwa_input_map["reference_fasta"]    : "ref{}.fa".format(uid1),
        bwa_input_map["reads_master_fasta"] : "{}.fa".format(uid2),
    }

    for suffix in bwa_index_file_suffixes():
        user_paths[bwa_index_map[suffix]] = "ref{uid}.fa{suff}".format(uid=uid1, suff=suffix)

    localFiles = LocalFileManager(job=job,
                                  fileIds_to_get=bwa_input_map.values() + bwa_index_map.values(),
                                  userFileNames=user_paths)

    dkr_reference_path = DOCKER_DIR + localFiles.localFileName(bwa_input_map["reference_fasta"])
    dkr_reads_path     = DOCKER_DIR + localFiles.localFileName(bwa_input_map["reads_master_fasta"])
    bwa_mem_parameters = ["mem", "-x", "ont2d", dkr_reference_path, dkr_reads_path]
    output_path        = localFiles.workDir() + "aln{}.sam".format(uid3)

    with open(output_path, 'w') as out_aln:
        tlp.docker_call(tool=bwa_docker_image,
                        parameters=bwa_mem_parameters,
                        work_dir=localFiles.workDir(),
                        outfile=out_aln)

    assert os.path.exists(output_path)
    alignment_fid = job.fileStore.writeGlobalFile(output_path)

    return AlignmentStruct(alignment_fid, AlignmentFormat.SAM)


def bwa_export_alignment(job, bwa_output_map, out_sam_path):
    job.fileStore.exportFile(bwa_output_map["alignment"], out_sam_path)
    return


def bwa_docker_alignment_root(job, config,
                              memory="10M", cores=1, disk="10M",
                              bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    # maps the various files needed to their unique fileStoreId, used
    # throughout the alignment pipeline
    # TODO this extra dict is unnecessary get rid of it
    bwa_input_map = {
        "reference_fasta"   : config["reference_FileStoreID"],
        "reads_master_fasta": config["sample_FileStoreID"],
    }
    bwa_index_map = job.addChildJobFn(bwa_index_docker_call, bwa_input_map).rv()
    alignment_job = job.addFollowOnJobFn(bwa_docker_align, bwa_input_map, bwa_index_map)
    bwa_alnstruct = alignment_job.rv()

    return bwa_alnstruct
