# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:28:39 2023

@author: chatt
"""
import numpy as np
import sys, re 
import math
from numpy import mean, sqrt, array 

def calc_RadGy(posArr):
    # posArr = N,3 array for N sites
    com = mean(posArr, axis=0) # center of mass
    Rg2 = mean(np.sum((posArr - com)**2, axis=1))
    return com, sqrt(Rg2)  

def parseXYZfile(file):
    with open(file, 'r') as tf:
        lines = tf.readlines()
        
    tmp_arr = array([line.split() for line in lines[2:]])
    return tmp_arr[:,1:].astype(np.float32)

def fixAtomNumbers(n, NA, NB, seg):
    atomNumbers = "			atomNumbers {}\n"
    chainLen = len(seg)*int(n)
    numAtoms = chainLen * (int(NA) + int(NB)) 
    for i,n in enumerate(range(1, numAtoms+1, chainLen//2)):
        if i == 0: continue
        if i == 1: 
            atomNumbers = atomNumbers[:-2] + str(n) + atomNumbers[-2:]
            continue
        atomNumbers = atomNumbers[:-2] + " " + str(n) + atomNumbers[-2:]
    return atomNumbers

def getBoxDim(posArr):
    x,y,z = posArr[:,0], posArr[:,1], posArr[:,2]
    return round(max(x)-min(x)), round(max(y)-min(y)), round(max(z)-min(z))

def updateCV(cvfile, Rg, L, n, NA, NB, seg):
    Rg = int(round(Rg))

    with open(cvfile, 'r') as tf:
        lines = tf.readlines()
    for i, line in enumerate(lines):
        if re.search('upperBoundary', line):
            lines[i] = f'\tupperBoundary {Rg + 10}\n'
        if re.search ('upperWalls', line):
            lines[i] = f'upperWalls {Rg + 5}\n'
        if re.search('atomNumbers', line):
            lines[i] = fixAtomNumbers(n, NA, NB, seg)
    
    #ofile = cvfile.split('.')[0]
    #LN = int(float(L)+2)
    ofile = f'N200_Rg_L{L}.colvars'

    with open(ofile, 'w') as tmpf:
        tmpf.writelines(lines)

_, xyzFile, cvfile, L, n, NA, NB, seg = sys.argv


#file = '//wsl.localhost/Ubuntu/home/achattaraj/simData_WSL/lammps_sims/lammps_sim_setup/IC_N160.xyz'

posArr = parseXYZfile(xyzFile)

com, Rg = calc_RadGy(posArr)

updateCV(cvfile, Rg, L, n, NA, NB, seg)

print()
print('File: ', xyzFile)

print(f'RadGy ~ {int(round(Rg))} A')

print('Dimension length (x,y,z): ', getBoxDim(posArr))


