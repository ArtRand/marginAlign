"""hmm.py Hidden Markov Model class
"""
import math
import uuid
import random

from toil_lib import require

from localFileManager import LocalFile, urlDownload

SYMBOL_NUMBER = 4


class Hmm:
    def __init__(self, modelType):
        self.modelType = modelType
        self.stateNumber = { "fiveState":5, "fiveStateAsymmetric":5, "threeState":3, "threeStateAsymmetric":3}[modelType]
        self.transitions = [0.0] * self.stateNumber**2
        self.emissions = [0.0] * (SYMBOL_NUMBER**2 * self.stateNumber)
        self.likelihood = 0.0
        self.runningLikelihoods = []

    def _modelTypeInt(self):
        return { "fiveState":0, "fiveStateAsymmetric":1, "threeState":2, "threeStateAsymmetric":3}[self.modelType]

    def write(self, file):
        f = open(file, 'w')
        f.write(("%s " % self._modelTypeInt()) + " ".join(map(str, self.transitions)) + (" %s\n" % self.likelihood))
        f.write(" ".join(map(str, self.emissions)) + "\n")
        f.close()

    def addExpectationsFile(self, file):
        fH = open(file, 'r')
        l = map(float, fH.readline().split())
        assert len(l) == len(self.transitions)+2
        assert int(l[0]) == self._modelTypeInt()
        self.likelihood += l[-1]
        self.transitions = map(lambda x : sum(x), zip(self.transitions, l[1:-1]))
        assert len(self.transitions) == self.stateNumber**2
        l = map(float, fH.readline().split())
        assert len(l) == len(self.emissions)
        self.emissions = map(lambda x : sum(x), zip(self.emissions, l))
        assert len(self.emissions) == self.stateNumber * SYMBOL_NUMBER**2
        self.runningLikelihoods = map(float, fH.readline().split()) #This allows us to keep track of running likelihoods
        fH.close()
        return self

    @staticmethod
    def loadHmm(file):
        fH = open(file, 'r')
        l = fH.readline().split()
        assert len(l) > 0
        fH.close()
        return Hmm({ 0:"fiveState", 1:"fiveStateAsymmetric", 2:"threeState", 3:"threeStateAsymmetric"}[int(l[0])]).addExpectationsFile(file)

    def normalise(self):
        """Normalises the EM probs.
        """
        for fromState in xrange(self.stateNumber):
            i = self.stateNumber * fromState
            j = sum(self.transitions[i:i+self.stateNumber])
            for toState in xrange(self.stateNumber):
                self.transitions[i + toState] = self.transitions[i + toState] / j
        for state in xrange(self.stateNumber):
            i = state * SYMBOL_NUMBER * SYMBOL_NUMBER
            j = sum(self.emissions[i:i+SYMBOL_NUMBER * SYMBOL_NUMBER])
            for emission in xrange(SYMBOL_NUMBER * SYMBOL_NUMBER):
                self.emissions[i + emission] = self.emissions[i + emission] / j

    def randomise(self):
        """Randomise the values in the HMM to small values.
        """
        self.transitions = map(lambda x : random.random(), range(self.stateNumber*self.stateNumber))
        self.emissions = map(lambda x : random.random(), range(self.stateNumber*SYMBOL_NUMBER*SYMBOL_NUMBER))
        self.normalise()

    def equalise(self):
        """Initialise the hmm with all equal probabilities.
        """
        self.transitions = [1.0/self.stateNumber] * (self.stateNumber**2)
        self.emissions = [1.0/(SYMBOL_NUMBER*SYMBOL_NUMBER)] * (self.stateNumber * SYMBOL_NUMBER**2)

    def setEmissionsToJukesCantor(self, divergence):
        i = (0.25 + 0.75*math.exp(-4.0*divergence/3.0))/4.0
        j = (0.25 - 0.25*math.exp(-4.0*divergence/3.0))/4.0
        for state in xrange(self.stateNumber):
            for x in xrange(SYMBOL_NUMBER):
                for y in xrange(SYMBOL_NUMBER):
                    self.emissions[(state * SYMBOL_NUMBER**2) + x * SYMBOL_NUMBER + y] = i if x == y else j

    def tieEmissions(self):
        """Sets the emissions to reflect overall divergence, but not to distinguish between different base identity
        """
        for state in xrange(self.stateNumber):
            a = self.emissions[state*SYMBOL_NUMBER**2:(state+1)*SYMBOL_NUMBER**2]
            identityExpectation = sum(map(lambda i : float(a[i]) if (i % SYMBOL_NUMBER) == (i / SYMBOL_NUMBER) else 0.0, xrange(SYMBOL_NUMBER**2)))
            a = map(lambda i : identityExpectation/SYMBOL_NUMBER if (i % SYMBOL_NUMBER) == (i / SYMBOL_NUMBER) else (1.0 - identityExpectation)/(SYMBOL_NUMBER**2 - SYMBOL_NUMBER), xrange(SYMBOL_NUMBER**2))
            assert sum(a) + 0.001 > 1.0 and sum(a) - 0.001 < 1.0
            self.emissions[state*SYMBOL_NUMBER**2:(state+1)*SYMBOL_NUMBER**2] = a
        assert len(self.emissions) == self.stateNumber * SYMBOL_NUMBER**2

    @staticmethod
    def modelFilename(global_config, get_url=False):
        model_filename = "{lab}_{emiter}_trainedmodel.hmm".format(lab=global_config["sample_label"],
                                                                  emiter=global_config["em_iterations"])
        if get_url:
            return global_config["output_dir"] + model_filename
        else:
            return model_filename


def downloadHmm(parent_job, config):
    if config["EM"]:
        hmm_file = LocalFile(workdir=parent_job.fileStore.getLocalTempDir(), filename=uuid.uuid4().hex)
        urlDownload(parent_job, Hmm.modelFilename(config, True), hmm_file)
        return Hmm.loadHmm(hmm_file.fullpathGetter())
    else:
        require(config["input_hmm_FileStoreID"],
                "[realignSamFileJobFunction]Need to provide a HMM for alignment or perform alignment/EM")
        hmm_file = parent_job.fileStore.readGlobalFile(config["input_hmm_FileStoreID"])
        return Hmm.loadHmm(hmm_file)
