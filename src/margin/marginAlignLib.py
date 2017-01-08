from margin.utils import *
import numpy as np


def normaliseHmmByReferenceGCContent(hmm, gcContent):
    """Normalise background emission frequencies to GC% given
    """
    for state in range(hmm.stateNumber):
        if state not in (2, 4): # Don't normalise GC content of insert states 
                                # (as they don't have any ref bases!)
            n = toMatrix(hmm.emissions[(SYMBOL_NUMBER**2) * 
                                       state:(SYMBOL_NUMBER**2) * (state+1)])
            hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = \
            fromMatrix(map(lambda i : map(lambda j : (n[i][j]/sum(n[i])) * 
            (gcContent/2.0 if i in [1, 2] else (1.0-gcContent)/2.0), range(SYMBOL_NUMBER)), 
                           range(SYMBOL_NUMBER))) #Normalise

def setHmmIndelEmissionsToBeFlat(hmm):
    """Set indel emissions to all be flat
    """
    for state in range(1, hmm.stateNumber):
        hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = \
        [1.0/(SYMBOL_NUMBER**2)]*SYMBOL_NUMBER**2  

def modifyHmmEmissionsByExpectedVariationRate(hmm, substitutionRate):
    #Normalise background emission frequencies, if requested to GC% given
    n = toMatrix(map(lambda i : (1.0-substitutionRate) if i % SYMBOL_NUMBER == \
                     i / SYMBOL_NUMBER else substitutionRate/(SYMBOL_NUMBER-1), 
                     xrange(SYMBOL_NUMBER**2)))
    hmm.emissions[:SYMBOL_NUMBER**2] = fromMatrix(np.dot(toMatrix(hmm.emissions[:SYMBOL_NUMBER**2]), n))

