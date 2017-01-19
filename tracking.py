
from sklearn.base import BaseEstimator
from sklearn.cluster import DBSCAN

import transform


class ClusterDBSCAN(BaseEstimator):
    def __init__(self):
        self.eps = 0.01
        self.min_hits = 5
        self.cls = DBSCAN(eps=self.eps, min_samples=self.min_hits)
    
    def fit(self, X, y):
        X=transform.polar(X)
        #        print X
        self.cls.fit(X,y)

    def predict(self, X):
        X=transform.polar(X)
        #        print X
        return self.cls.fit_predict(X)

