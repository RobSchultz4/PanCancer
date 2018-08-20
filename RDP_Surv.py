#################################################################################################################
###
### Replicate survival
###
#################################################################################################################

#This will be the survival code implemented into main when it is working and generalized. for now, I will try to get it to work on one case.

from lifelines import CoxPHFitter
import csv
import pandas as pd
import numpy as np

def runcph(coxdf, dur, eve):
    cph = CoxPHFitter()
    cph.fit(coxdf, duration_col = dur, event_col = eve, show_progress=True)  
    return cph

def isincluded(p1, colname):
    return any(pd.isnull(p1[colname]).apply(lambda x: not x))

def areValues(p1, colname):
    return  pd.isnull(p1[colname]).apply(lambda x: not x)

def makeMask(colname1,colname2):
    mask=[True for i in  np.arange(0,len(p1[colname1]))]
    mask1 = areValues(p1,colname1)
    mask2 = areValues(p1,colname2)
    for i in np.arange(0,len(p1[colname1])):
        if mask1[i] == False:
            mask[i] = False
        elif mask2[i] == False:
            mask[i] = False
        else:
            mask[i] = True 
    return mask

def handlenans(p1, colname1, colname2, coxdf):
    if any(pd.isnull(p1[colname1])) or any(pd.isnull(p1[colname2])): 
        mask = makeMask(colname1,colname2)
        coxdf = coxdf[mask]
        status = p1[colname2][mask]
        time = p1[colname1][mask]
    else:
        time = p1[colname1]
        status = p1[colname2]
    return time, status, coxdf

p1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed pData/GSE49278_pData_recon.csv',header = 0,index_col = 0)
coxdfo = pd.DataFrame([])
coxdfp = pd.DataFrame([])
# Run if there is any survival info availible
if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME') or isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
    coxdfo['Bic Expression'] = pcbic
    coxdfp['Bic Expression'] = pcbic
    # Run if the info availible is Overall Survival
    if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME'):
        time, status, coxdfo = handlenans(p1,'OS.STATUS','OS.TIME', coxdfo)
        coxdfo['OS.STATUS'] = status
        coxdfo['OS.TIME'] = time
        survo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
        # Run if there is also age information
        if isincluded(p1, 'AGE'):
            coxdfo['AGE'] = p1['AGE']
            survAo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
            # Run if there is also sex information
            if isincluded(p1, 'SEX'):
                coxdfo['SEX'] = pd.get_dummies(p1['SEX'])['M'] # Male is 1 and female is 0
                survASo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
    # Run if there is Progression Free Survival
    if isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
        time, status, coxdfp = handlenans(p1,'PFS.STATUS','PFS.TIME', coxdfp)
        coxdfp['PFS.STATUS'] = status
        coxdfp['PFS.TIME'] = time 
        survp=runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
        # Run if there is also age information
        if isincluded(p1, 'AGE'):
            coxdfp['AGE'] = p1['AGE']
            survAp=runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
            # Run if there is also sex information
            if isincluded(p1, 'SEX'):
                coxdfp['SEX'] = pd.get_dummies(p1['SEX'])['M'] # Male is 1 and female is 0
                survASp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
else:
    print("No Survival!")
    






