#!/usr/bin/env python2.7
# oogenerateParticles.py
# Thomas Boser

"""
Note to self: particle printing format
particle barcode, [ vertex x,y,z], [ momentum p, theta, phi ], charge
4509578221846528, [-0.0224014, -0.00570905, -51.0684], [1450.18, 2.06558, -2.19591], 1
"""

from __future__ import print_function, division

import random
import math
import oogenerateHits as gh
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sympy.geometry import Circle, Ray, Point, intersection 

class Particle:
    """ particle constructor """
    def __init__(self, barcode, gen = False, vertices=[]):
        self.barcode = barcode
        self.vertices = vertices
        self.mangle = []
        self.charge = -10

        self.hits = []
        self.hitbcs = []

        ### FOR GENERATED PARTICLES ###
        self.p_radius = 0
        self.centpt = []
        self.p_circ = None
        self.m_ray = None

        if(gen): self.genParticle()

    def genParticle(self):
        """ generate values for particle """
        self.vertices = [round(random.uniform(.03, -.03),7), round(random.uniform(.03, -.03),7)\
                         ,0] #Leaving z at 0, only considering 2d space for now
        self.mangle = [round(random.uniform(15000, 4000),2), round(random.uniform(4, 0),5),\
                       round(random.uniform(4, -4),5)]
        self.charge = random.randint(-1, 1)

        #particle trajectory estimations
        magfield = 1 # we will assume the magnetic field is uniform
        if self.charge != 0: 
            self.p_radius = self.mangle[0] / (self.charge * magfield)
            self.centpt = [self.vertices[0] + (self.p_radius*math.sin(self.mangle[2])),
                           self.vertices[1] + (self.p_radius*math.cos(self.mangle[2]))]
            self.p_circ = Circle(Point(self.centpt[0], self.centpt[1]), abs(self.p_radius))
        self.m_ray = Ray(Point(self.vertices[0], self.vertices[1]), #ray headed in direction of momentum (angle phi)
                           Point(100, 100*math.tan(self.mangle[2])))

    def addHit(self, hit):
        """ add a hit to self.hits """
        if hit not in self.hits: self.hits.append(hit)

    def getHits(self, detectors, hitbc):
        """ produces all hits with the detectors specified, appends them
            to self.hits, assumes particles travel in circular motion so
            this is an approximation """
        if len(self.hits) == 0: 
            for detpos, det in enumerate(detectors):
                poss_hits = self.getIntersects(det)
                if poss_hits is not None:
                    hitbc += 1
                        #                  print("creating hit: ",hitbc, self.barcode, poss_hits,detpos+1)
                    self.hits.append(gh.Hit(hitbc, self.barcode, poss_hits,detpos+1))
        return hitbc

    def getIntersects(self, detector):
        """ returns intersection pts of two cirles (or a line and a circle).
            does so using the sympy library. """
        m_intersect = detector.circle.intersection(self.m_ray)
        if self.charge == 0: #if charge is 0 particle travels in straight line
            return np.array((m_intersect[0].x.evalf(), m_intersect[0].y.evalf()))
        else: #circle circle intersection if particle is charged
            intersects = detector.circle.intersection(self.p_circ)
            if len(intersects) == 0: return None
            if len(intersects) == 1: 
                return np.array((intersects[0].x.evalf(), intersects[0].y.evalf()))
            flag = 0
            if intersects[0].distance(m_intersect[0]) > intersects[1].distance(m_intersect[0]):
                flag = 1
            return np.array((intersects[flag].x.evalf(), intersects[flag].y.evalf()))

    def getAngle(self):
        """ return the angle of the particle's trail """
        if self.charge > 0: return 0.0
        return 180.0

    def getTheta1(self):
        """ return theta1 for the particle arc """
        angle = 0
        if self.vertices[0] - self.hits[-1].lhit[0] < 0 or self.vertices[1] - self.centpt[1] < 0: 
            angle = math.degrees(math.atan((self.hits[-1].lhit[1] - self.centpt[1]) /\
                                           (self.hits[-1].lhit[0] - self.centpt[0]))) #flipped thetas bc arc drawn counter clockwise
            if self.centpt[0]- self.hits[-1].lhit[0] > 0:
                angle = 180 + abs(angle)
        else: 
            angle = math.degrees(math.atan((self.vertices[1] - self.centpt[1]) /\
                                           (self.vertices[0] - self.centpt[0])))
            if self.centpt[0] - self.vertices[0] > 0:
                angle = 180 - abs(angle)
        return angle

    def getTheta2(self):
        """ return theta2 for the particle arc """
        angle = 0
        if self.vertices[0] - self.hits[-1].lhit[0] < 0 or self.vertices[1] - self.centpt[1] < 0:
            angle = math.degrees(math.atan((self.vertices[1] - self.centpt[1]) /\
                                           (self.vertices[0] - self.centpt[0])))
            if self.centpt[0] - self.vertices[0] > 0:
                angle = 180 + abs(angle)
        else: 
            angle = math.degrees(math.atan((self.hits[-1].lhit[1] - self.centpt[1]) /\
                                           (self.hits[-1].lhit[0] - self.centpt[0])))
            if self.centpt[0]- self.hits[-1].lhit[0] > 0:
                angle = 180 - abs(angle)
        return angle

    def plotOrigin(self):
        """ plot the origin point for the particle """
        plt.scatter(self.vertices[0], self.vertices[1], c='r')

    def plotTrail(self):
        """ plot the trail the particle follows """ 
        if len(self.hits) == 0: return #if there are no hits there is no trail
        if not self.charge == 0: #plot an arc
            ax = plt.gca()
            print("angle =", self.getAngle(), "theta1 =", 
                  self.getTheta1(), "theta2 =", self.getTheta2()) ###TEMP
            ax.add_patch(patches.Arc(self.centpt, self.p_radius*2, self.p_radius*2,
                                     angle=0.0, theta1=self.getTheta1(),
                                     theta2=self.getTheta2()))
        else: #plot a line if charge is 0
            plt.plot((self.vertices[0], self.hits[-1].lhit[0]), 
                     (self.vertices[1], self.hits[-1].lhit[1]),
                      color = random.choice('bgrcmykw'))

    def plotJoins(self):
        """ plot all hits joined by a line """
        col = random.choice('bgrcmyk') #chose random color
        #self.hits.sort(key=lambda x: x.detpos)
        origin = self.vertices
        for hit in self.hits:
            plt.plot((origin[0], hit.lhit[0]),
                     (origin[1], hit.lhit[1]),
                      color = col)
            origin = hit.lhit

    def plotHits(self):
        """ plot all hits in self.hits """
        for hit in self.hits:
            hit.plotHit()

    def printParticle(self):
        """ print particle to stdout """
        print(self.barcode,',',self.vertices,',',self.mangle,',',self.charge, sep='')

    def printHits(self):
        """ print hits to stdout """
        for hit in self.hits:
            hit.printHit()

