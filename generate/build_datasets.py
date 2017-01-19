#!/usr/bin/env python2.7
# fileControl.py
# Thomas Boser

""" 
Creates datasets.
"""

import os
import sys
import particleController as pc

def createDataset(n):
    #first and second datasets
    orig_stdout = sys.stdout
    for file in ["hits_train", "particle_valid", "particle_test"]:
        cont = pc.particleController()
        cont.clearController()
        cont.createParticles(n)
        for i in range(1000, 10001, 1000): cont.addDetector(i)
        cont.computeallHits()
        cont.moveHits()
        outf = open(file + ".csv", 'w')
        sys.stdout = outf
        if file is "hits_train.csv": cont.printallHits()
        else: cont.printallHits(dataset = True)
        outf.close()
        outf = open(file + "_soln.txt", 'w')
        sys.stdout = outf
        cont.printSoln()
        outf.close()
    #zip input data, clean directory
    #    os.system("zip input_data_1_2 hits_train.csv particle_valid.csv particle_test.csv")
    #    os.system("rm hits_train.csv particle_valid.csv particle_test.csv")
    #    os.system("mkdir data_results")
    #    os.system("mv hits_train_soln.txt particle_valid_soln.txt particle_test_soln.txt data_results")

createDataset(100)
