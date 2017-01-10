#!/usr/bin/env python
import re
import os
import urlparse
import uuid
import subprocess
from toil_lib import require
from toil_lib.files import copy_files
from toil_lib.programs import docker_call


DEBUG = True
DOCKER_DIR = "/data/"

def removeTempSuffix(filename):
    # type: (string)
    """removes '.tmp' from a filename string
    """
    assert filename.endswith(".tmp")
    return re.sub('\.tmp$', '', filename)


class LocalFileManager(object):
    """Gets all the files in 'fileIds_to_get' from the global fileStore
    and puts them in the local working directory.
    returns: file_dict<file_id, (full_path, unique_file_name)>, work_dir
    """
    def __init__(self, job, fileIds_to_get, userFileNames=None):
        # type: (toil.job.Job, list<str>, dict<str, str>)
        # userFileNames has the FileStoreID as keys and file name  
        # you want it to have as a value {fid: "file_name"}
        self.owner_job = job
        self.work_dir  = job.fileStore.getLocalTempDir()
        self.file_dict = self._initialize_file_dict(fileIds_to_get, userFileNames)

    def localFileName(self, fileId):
        return self.file_dict[fileId][1]

    def localFilePath(self, fileId):
        return self.file_dict[fileId][0]

    def localFileUrl(self, fileId):
        return "file://" + self.file_dict[fileId][0]

    def workDir(self):
        return self.work_dir + "/"

    @staticmethod
    def safeLocalDelete(local_file_list):
        removed_bools = []
        paths = [urlparse.urlparse(f).path for f in local_file_list]
        for url in paths:
            assert os.path.exists(url), "Didn't find {}".format(url)
            os.remove(url)
            removed_bools.append(not os.path.exists(url))
        return removed_bools

    def _initialize_file_dict(self, file_ids, userFileNames):
        # TODO make this a comprehension
        file_dict = {}
        for f in file_ids:
            file_dict[f] = ""  # TODO make this None?

        for file_id in file_dict.keys():
            if DEBUG:
                self.owner_job.fileStore.logToMaster("[LocalFileManager]Getting {}".format(file_id))
            if userFileNames is None:
                temp_file_name = uuid.uuid4().hex + ".tmp"
            else:
                assert file_id in userFileNames.keys(), \
                    "Didn't find user-specified path for {}".format(file_id)
                temp_file_name = userFileNames[file_id]
            destination_path = self.workDir() + "{fn}".format(fn=temp_file_name)
            self.owner_job.fileStore.readGlobalFile(file_id, userPath=destination_path)
            assert(os.path.exists(destination_path)), "[LocalFileManager] Error getting file {}".format(file_id)
            file_dict[file_id] = (destination_path, temp_file_name)

        return file_dict


class LocalFile(object):
    """A struct containing the path and handle for a file, used to easily access the a file and
    it's contents primarly useful for files that aren't in the FileStore already
    """
    def __init__(self, workdir, filename):
        # TODO make this an os.path.join
        self.path     = workdir + "/" + filename
        self.filename = filename
        self.workdir  = workdir

    def filenameGetter(self):
        return self.filename

    def fullpathGetter(self):
        return self.path

    def workdirGetter(self):
        return self.workdir


def deliverOutput(parent_job, deliverable_file, destination, retry_count=3, s3am_image="quay.io/ucsc_cgl/s3am"):
    # type (toil.job.Job, LocalFile, URL)
    """Run S3AM in a container to deliver files to S3 or copy files to a local file
    """
    require(os.path.exists(deliverable_file.fullpathGetter()), "[deliverOutput]Didn't find file {here}"
                                                               "".format(here=deliverable_file))

    if urlparse.urlparse(destination).scheme == "s3":
        parent_job.fileStore.logToMaster("[deliverOutput]Using S3AM docker...")
        source_arg = DOCKER_DIR + deliverable_file.filenameGetter()
        s3am_args  = ["upload", source_arg, destination]
        for i in xrange(retry_count):
            try:
                docker_call(tool=s3am_image, parameters=s3am_args, work_dir=(deliverable_file.workdirGetter() + "/"))
            except subprocess.CalledProcessError:
                parent_job.fileStore.logToMaster("[deliverOutput]S3AM failed with args {}".format(s3am_args.__str__()))
            else:
                parent_job.fileStore.logToMaster("[deliverOutput]S3AM succeeded")
                return
        raise RuntimeError("[deliverOutput]Delivering {deliverable} to {destination} failed after {n} attempts"
                           "".format(deliverable=deliverable_file, destination=destination, n=retry_count))
    else:
        require(urlparse.urlparse(destination).scheme == "file", "[deliverOutput]Illegal output URL {}"
                                                         "".format(destination))
        destination_dir = urlparse.urlparse(destination).path
        parent_job.fileStore.logToMaster("[deliverOutput]copying sam {f} to output {out}"
                                         "".format(f=deliverable_file.fullpathGetter(), out=destination_dir))
        copy_files(file_paths=[deliverable_file.fullpathGetter()], output_dir=destination_dir)
        return
