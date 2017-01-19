import pandas as pd
import numpy as np

from sklearn.cross_validation import ShuffleSplit

import tracking
import score_cluster
import score_physics
import plotting

filename = "test.csv"

def read_data(filename):
    df = pd.read_csv(filename)
    y_df = df[['particle']]
    X_df = df.drop(['hit','particle'], axis=1)
    return X_df, y_df



if __name__ == '__main__':
    print("Reading file ...")

    X_df, y_df = read_data(filename)

    skf = ShuffleSplit(
    len(y_df), n_iter=1, test_size=0.2, random_state=57)
    print("Training file ...")
    for train_is, test_is in skf:
        print '--------------------------'
#        print train_is
#        print test_is

        tracker = tracking.ClusterDBSCAN()

        X_train_df = X_df.iloc[train_is].copy()
        y_train_df = y_df.iloc[train_is].copy()
        X_test_df = X_df.iloc[test_is].copy()
        y_test_df = y_df.iloc[test_is].copy()

        # Temporarily bypass splitting (need to avoid shuffling events)
        X_train_df = X_df.copy()
        y_train_df = y_df.copy()
        X_test_df = X_df.copy()
        y_test_df = y_df.copy()
        
#        print X_train_df
        tracker.fit(X_train_df.values, y_train_df.values)

        y_predicted = tracker.predict(X_test_df.values)
        score_eff, score_fake = score_cluster.evaluate(y_test_df.values, y_predicted)
        print 'score efficiency = ', score_eff
        print 'score fake = ', score_fake
