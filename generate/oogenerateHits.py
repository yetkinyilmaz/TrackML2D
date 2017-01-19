#!/usr/bin/env python2.7
# oogenerateHits.py
# Thomas Boser

"""
Note to self: hit printing format
hit barcode, particle barcode, [ local hit x,y ] , [ local error ex, ey ], [ global x, y, z ], [ [ channelX, channlY, chargeInChannel ]]
578, 27022972153757696, [ -32.9653, -26.4], [ 0, 0],[ -167.124, -393.564,1553.25], [[ 139, 13, 0.0249668], [ 140, 13, 0.286339], [ 141, 13, 0.105968]]
"""

from __future__ import print_function

import math
import numpy as np
import matplotlib.pyplot as plt

class Hit:
    """ hit constructor """
    def __init__(self, hbc, pbc, point, detpos=0):
        self.hbc = int(hbc)
        self.pbc = pbc
        self.lhit = point
        self.detpos = detpos #to speed up certain processes

        ### CURRENTLY UNUSED ###
        self.lerr = []
        self.globhit = []
        self.cxcycharge = [[]*3]

    def __str__(self):
        """ overload print operator """
        return str(self.hbc)

    def __repr__(self):
        """ another print overload """
        return str(self.hbc)

    def origDist(self):
        """ return distance from origin """
        return math.sqrt(self.lhit[0]**2 + self.lhit[1]**2)

    def plotHit(self):
        """ plot a hit (single point) to plt figure """
        plt.scatter(self.lhit[0], self.lhit[1])

    def printHit(self, dataset = False):
        """ print hit to stdout """
        if dataset:
            print(self.hbc, ',', self.pbc, ',', list(self.lhit), sep='')
