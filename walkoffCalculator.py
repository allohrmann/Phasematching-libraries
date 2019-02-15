# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 14:52:34 2019

@author: cqt
"""

import nonlinearCrystal
import helper


pump = 0.405
s= 0.81

targetDisplacement = 1e3 # everything always in microns, except for sometimes randomly not
LBBO1 = 10e3 # initial guess
LBBO2 = 11e3 #initial guess
opticalAxis = 45 # in degrees
displacementTolerance = 0.1 # not really achievable, but ok

BBOI = nonlinearCrystal.NonlinearCrystal(l = LBBO1, material = "BBO", theta =opticalAxis)
BBOII = nonlinearCrystal.NonlinearCrystal(l = LBBO2, material = "BBO", theta =opticalAxis)



def gradientDescent(target,crystal, wl):
    
    slope = 1
    d1 = crystal.calcWalkoff(wl, opticalAxis)*crystal.l
    #d2 = BBOII.calcWalkoff(0.81, opticalAxis)*BBOII.l

    grad = d1 - target
    
    while abs(d1 - target) >  displacementTolerance:
        
        grad = (d1 - target) * slope
        crystal.l -= grad
        d1 = crystal.calcWalkoff(wl, opticalAxis)*crystal.l



gradientDescent(targetDisplacement,BBOI, pump)
gradientDescent(targetDisplacement, BBOII, s)

print ("###########################################")
d1 = (BBOI.calcWalkoff(0.405, opticalAxis)*BBOI.l)
d2 = (BBOII.calcWalkoff(0.81, opticalAxis)*BBOII.l)

print ("BBO1 legnth:   " + str(BBOI.l/1000) + " mm")
print ("BBO2 length:   " + str(BBOII.l/1000) + " mm")
            

print ("Walkoff splitter:   " + str(d1))
print ("Walkoff combiner:   " + str(d2))
            
print ("Walkoff difference: " + str(d1-d2))
print ("Relative walkoff difference:  " + str((d1-d2)/d1) )
print ("###########################################")



