# -*- coding: utf-8 -*-

import os
import numpy as np
import math
import IPython

solname = str(input('Enter solvent name: \n'))
K = 3

modelID = []
for line in open(os.path.join('keys', 'models.txt')):
    if line.startswith('D'):
        info=line.split()
        modelID.append(info[0])

AICdict = {}
for k in range(len(modelID)):
    delGname = "delG_HHR-%s.dat" % modelID[k]
    yintname = "yint-%s-%s.dat" % (modelID[k], solname)
    
    calHHR = []
    exphyd = []
    calhyd = []
    
    for line in open(os.path.join('data', modelID[k], delGname)):
        #Might need to fix this part later if we end up having more than 999 "systems"
        if line.startswith('Y0'):
            info=line.split()
            calHHR.append(float(info[5]))
            exphyd.append(float(info[6]))
            
    yintfile = open(os.path.join('data', modelID[k], yintname), "r")
    yint = float(list(yintfile.readlines())[0].split()[0])
    yintfile.close()
    
    unc = 0
    for i in range(len(calHHR)):
        calhyd.append(calHHR[i]-yint)
        unc += (exphyd[i] - calhyd[i])**2
    unc = unc/len(calHHR)
    print("Model ID: ", modelID[k], "\nOriginal best estimate for uncertainty: ", unc)
    
    unclist = []
    for fluct in np.linspace(-1,1,50):
        unclist.append(unc+fluct*unc/16)
    lnLlist = []
    for i in range(len(unclist)):
        lnL = 0
        for j in range(len(calHHR)):
            lnL += (exphyd[j] - calhyd[j])**2 / unclist[i] + 2*np.log(math.sqrt(2*math.pi*unclist[i]))
        lnL = -1/2 * lnL
        lnLlist.append(lnL)
    maxlnL = max(lnLlist)
    maxindex = lnLlist.index(maxlnL)
    
    unc = unclist[maxindex]
    print("Final uncertainty: ", unc, "\nIndex # (out of %s): " % len(unclist), maxindex, "\n")
    
    AIC = -2*maxlnL + 2*K + (2*K*(K+1))/(len(calHHR)-K-1)
    AICdict[modelID[k]] = AIC
    
minAICkey = min(AICdict, key=AICdict.get)
relprobdict = {}
for k in range(len(modelID)):
    relprob = math.e ** ((AICdict[minAICkey]-AICdict[modelID[k]])/2)
    relprobdict[modelID[k]] = relprob

with open(os.path.join("data", "AIC_%s.dat" % solname), "w") as f:
    f.write("MID       AIC                     exp((AIC_min - AIC_i)/2)\n")
    f.write("---       ---                     ------------------------\n")
    for i in range(len(modelID)):
        f.write('{0:<10}{1:<24}{2:<24}\n'.format(modelID[i], AICdict[modelID[i]], relprobdict[modelID[i]]))
            
