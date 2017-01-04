class AlignmentFormat:
    SAM, BAM = range(2)


class AlignmentStruct(object):
    def __init__(self, filestoreId, alignment_format):
        # type (str, AlignmentFormat)
        self.fid        = filestoreId
        self.aln_format = alignment_format

    def FileStoreID(self):
        return self.fid

    def FID(self):
        return self.fid

    def AlignmentFormat(self):
        return self.aln_format

    def IsSAM(self):
        return (self.aln_format == AlignmentStruct.SAM)

    def IsBAM(self):
        return (self.aln_format == AlignmentFormat.BAM)
