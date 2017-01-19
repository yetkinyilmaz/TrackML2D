import numpy as np

def polar(hits, rscale = 0.0000001):
    ptsnew = np.zeros(hits.shape)
    xy = hits[:,0]**2 + hits[:,1]**2
    ptsnew[:,0] = np.sqrt(xy)*rscale
    ptsnew[:,1] = np.arctan2(hits[:,1], hits[:,0])
    return ptsnew


