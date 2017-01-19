import numpy as np
# import pandas as pd


def evaluate(y_test, y_pred):

    eff = 0;
    fake = 0.;
    
# evaluate noise
    for track in range(0,len(np.unique(y_pred))):
        fake+=0


# remove combinatorials


# evaluate hit efficiency
    for particle in np.unique(y_test):
        print "particle : ", particle
        true_hits = y_test[y_test[:,0] == particle]
        #        print "true hits : ", true_hits
        found_hits = y_pred[y_test[:,0] == particle]
        #        print "found hits : ", found_hits

        nsubcluster=len(np.unique(found_hits[found_hits[:] >= 0]))
        #        print "found clusters : ", nsubcluster

        if(nsubcluster > 0):
            maxcluster = np.argmax(np.bincount(found_hits[found_hits[:] >= 0]))
            eff = len(found_hits[found_hits[:] == maxcluster])/len(true_hits)

        print "efficiency : ", eff



    return eff, fake



