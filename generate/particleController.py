#!/usr/bin/env python2.7
# particleController.py
# Thomas Boser

from __future__ import print_function, division

import sys
import ast
import oogenerateParticles as gp
import oogenerateHits as gh
import oogenerateDetectors as gd
import matplotlib.pyplot as plt
from math import sqrt
from sympy.geometry import Circle, Ray, Point, intersection 

class particleController:
    """ controller class for particles, hits, and detectors """
    def __init__(self):
        """ controller constructor """
        self.bcTracker = 0 #track barcodes. doesn't matter for our purposes just can't have two of the same
        self.hbcTracker = 0

        self.particles = []
        self.detectors = []
        self.hits = []

        #weird way to track plot state, i'm sure there is a better way to do this
        self.plotState = []

    ### PARTICLE & HIT METHODS ###
    #-- PARTICLE CONTROL METHODS --#
    def addParticle(self, *args):
        """ add a particle to self.particles with preset parameters """
        for arg in args: self.particles.append(arg)

    def clearParticles(self):
        """ delete every particle in self.particles """
        self.particles[:] = []

    #-- HIT CONTROL METHODS --#
    def addHit(self, hbc, point, pbc="?"):
        """ read hits from files """
        self.hits.append(gh.Hit(hbc, pbc, point))

    def clearHits(self):
        self.hits[:] = []

    #-- PARTICLE CREATION METHODS --#
    def createParticles(self, n):
        """ add n particles to self.particles, increment bcTracker """
        for i in range(n):
            self.bcTracker += 1
            self.particles.append(gp.Particle(self.bcTracker, gen=True))

    def computeallHits(self):
        """ compute hits for all particles """
        for particle in self.particles:
            self.hbcTracker = particle.getHits(self.detectors, 
                                               self.hbcTracker) #update barcode tracker

    #-- IO METHODS --#
    def readHits(self, *args):
        """ read particles from files """
        for infilen in args:
            try:
                infile = open(infilen, 'r')
            except IOError:
                print("infile could not be opened")
                exit(1)
            for line in infile:
                line = line.replace(']','').replace('[','')\
                    .replace(' ','').rstrip().split(',')
                self.addHit(line[0], [float(line[2]), float(line[3])], line[1])

    def printParticles(self):
        """ print all particles in self.particles to stdout """
        for partc in self.particles:
            partc.printParticle()

    def printnumParticles(self):
        """ print number of particles initialized """
        print("There are", len(self.particles), "particles.")

    def printallHits(self, dataset = False):
        """ print all hits in self.hits to stdout """
        if len(self.hits) == 0: self.moveHits()
        for hit in self.hits:
            hit.printHit(dataset)

    def printnumHits(self):
        """ print number of computed hits """
        print("There are", len(self.hits), "hits.")

    def writeParticles(self):
        """ write all particles to a file """
        pass

    def printSoln(self):
        """ print solution to stdout """
        for particle in self.particles:
            print(particle.barcode, ",", particle.hits)

    #-- PARTICLE PREDICTION METHODS --#
    def predictParticles(self):
        """ predict particles using hits in self.hits """
        temphits = self.hits
        for hit in self.hits:
            if hit.pbc != "?":
                #note: assuming particle origin is (0, 0, 0)
                thispart = gp.Particle(hit.pbc, False, [0, 0])
                self.particles.append(thispart) # create particle
                thispart.hits.append(hit)
                thispart.hitbcs.append(hit.hbc)
                temphits.remove(hit)
        for particle in self.particles:
            for i in range(1, len(self.detectors)):
                minlen = 99999
                minind = 0
                if len(particle.hits) == 1:
                    thisray = Ray(Point(0, 0), Point(particle.hits[0].lhit[0],
                                                     particle.hits[0].lhit[1]))
                else:
                    l = len(particle.hits)
                    thisray = Ray(Point(particle.hits[l-2].lhit[0], particle.hits[l-2].lhit[1]),
                                  Point(particle.hits[l-1].lhit[0], particle.hits[l-1].lhit[1]))
                intersect = self.detectors[i].circle.intersection(thisray)
                if len(intersect) == 0: continue 
                thisi = [intersect[0].x.evalf(), intersect[0].y.evalf()]
                for ind, hit in enumerate(temphits):
                    if hit.detpos == i+1:
                        thislen = self.pt_dist(thisi, hit.lhit)
                        if thislen < minlen:
                            minlen = thislen
                            minind = ind
                particle.hits.append(temphits[minind])
                particle.hitbcs.append(temphits[minind].hbc)
                temphits[minind].pbc = particle.barcode
                #temphits.remove(temphits[minind])

    def compHitdet(self):
        """ specify which detector a hit intercepts """
        self.detectors.sort(key=lambda x: x.radius)
        for hit in self.hits:
            if hit.detpos == 0:
                for i, detector in enumerate(self.detectors):
                    if(self.isclose(detector.radius, hit.origDist(), 1)):
                        hit.detpos = i+1
                        break

    def scoreSoln(self, solnfile):
        """ score the result of self.predictParticles() """
        infile = open(solnfile, 'r')
        solutions = []
        for line in infile:
            solutions.append(ast.literal_eval(line))
        numcorrect = 0
        for solution in solutions:
            for particle in self.particles:
                if int(solution[0]) == int(particle.barcode):
                    numcorrect += len([val for val in solution[1] if int(val) in particle.hitbcs])
                    break
        print("You predicted ",(numcorrect/len(self.hits))*100, " '%' of hits correctly.", sep='')

    ### DETECTOR METHODS ###
    def addDetector(self, r):
        """ add detector of radius r to self.detectors """
        if gd.Detector(r) not in self.detectors: #prevent duplicate detectors
            self.detectors.append(gd.Detector(r))

    def clearDetectors(self):
        """delete every detector in self.detectors """
        self.detectors[:] = []

    def printDetectors(self):
        """ print all detectors in self.detectors to stdout """
        for det in self.detectors:
            det.printDetector()

    def printnumDetectors(self):
        """ print number of initialized detectors """
        print("There are", len(self.detectors), "detectors.")

    ### PLOTTING METHODS ###
    def plotParticle_trails(self):
        """ plot all particle trails """
        if self._plotParticle_trails not in self.plotState:
            self.plotState.append(self._plotParticle_trails)
        self.showPlot()

    def plotParticle_joins(self):
        """ plots line between points """
        if self._plotParticle_joins not in self.plotState:
            self.plotState.append(self._plotParticle_joins)
        self.showPlot()

    def plotParticle_origins(self):
        """ plot all particle origins """
        if self._plotParticle_origins not in self.plotState:
            self.plotState.append(self._plotParticle_origins)
        self.showPlot()

    def plotParticle_hits(self):
        """ plot all particle hits in self.particles.hits """
        if self._plotParticle_hits not in self.plotState:
            self.plotState.append(self._plotParticle_hits)
        self.showPlot()

    def plotallHits(self):
        """ plot all particle hits in self.hits """
        if self._plotallHits not in self.plotState:
            self.plotState.append(self._plotallHits)
        self.showPlot()

    def plotDetectors(self):
        """ plot all detectors """ 
        if self._plotDetectors not in self.plotState:
            self.plotState.append(self._plotDetectors)
        self.showPlot() 

    def clearPlot(self):
        """ clear plot """
        self.plotState[:] = []

    ### HELPER METHODS ###
    def moveHits(self):
        """ move hits from particles to particlecontroller after computing """
        self.clearHits()
        for particle in self.particles:
            self.hits.extend(particle.hits)

    def pt_dist(self, p1, p2):
        """ distance between two points described by a list """
        return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

    def isclose(self, v1, v2, tol):
        """ check if two values are within a tolerance range of eachother """
        return abs(v1 - v2) < tol

    def _plotParticle_origins(self):
        """ actually plots all particle origins, prevents recursion """
        for particle in self.particles:
            particle.plotOrigin()

    def _plotParticle_trails(self):
        """ actually plots all particle trails, prevents recursion """
        for particle in self.particles:
            particle.plotTrail()

    def _plotParticle_joins(self):
        """ actually plots line between points """
        for particle in self.particles:
            particle.plotJoins()

    def _plotallHits(self):
        """ actually plots all particle hits in self.hits, prevents recursion """
        for hit in self.hits:
            hit.plotHit()

    def _plotParticle_hits(self):
        """ actually plots all particle hits, prevents recursion """
        for particle in self.particles:
            particle.plotHits()

    def _plotDetectors(self):
        """ actually plots all detectors """
        for detector in self.detectors:
            detector.plotDetector()

    def showPlot(self):
        """ helper command to create new plot with commands used """
        for command in self.plotState: command()
        plt.scatter(0, 0, c='g') #always plot origin in green
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    ### CONTROLLER METHODS ###
    def clearController(self):
        """ clear controller of all values """
        self.clearDetectors()
        self.clearParticles()
        self.clearPlot()
        self.clearHits()