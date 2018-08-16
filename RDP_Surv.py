####################################################################################################################
###
### Replicate survival
###
####################################################################################################################

#This will be the survival code implemented into main when it is working and generalized for now, I will try to get it to work on one case.

from lifelines import CoxPHFitter
import csv
import pandas as pd
import numpy as np

def runcph(coxdf, dur, eve):
    coxdf = pd.DataFrame(coxdf).T
    cph = CoxPHFitter()
    cph.fit(coxdf, duration_col = dur, event_col = eve, show_progress=True)  
    return cph

def isincluded(p1, colname):
    return any(pd.isnull(p1[colname]).apply(lambda x: not x))

p1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed pData/GSE49278_pData_recon.csv',header = 0,index_col = 0)
coxdf = []

# Run if there is any survival info availible
if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME']) or isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
    # Run if the info availible is Overall Survival
    if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME']):
        coxdf.append(p1['OS.STATUS'])
        coxdf.append(p1['OS.TIME'])
        dur = 'OS.TIME'
        eve = 'OS.STATUS' 
        surv=runcph(coxdf, dur, eve)
        # Run if there is also age information
        if isincluded(p1, 'AGE'):
            coxdf.append(p1['AGE'])
            survA=runcph(coxdf, dur, eve)
            # Run if there is also sex information
            if isincluded(p1, 'SEX'):
                coxdf.append(pd.get_dummies(p1['SEX'])['M']) # Male is 1 and female is 0
                survAS=runcph(coxdf, dur, eve)
            
    # Run if there is Progression Free Survival
    if isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
        coxdf.append(p1['PFS.STATUS'])
        coxdf.append(p1['PFS.TIME'])
        dur = 'PFS.TIME'
        eve = 'PFS.STATUS'
        surv=runcph(coxdf, dur, eve)
        # Run if there is also age information
        if isincluded(p1, 'AGE'):
            coxdf.append(p1['AGE'])
            survA=runcph(coxdf, dur, eve)
            # Run if there is also sex information
            if isincluded(p1, 'SEX'):
                coxdf.append(pd.get_dummies(p1['SEX'])['M']) # Male is 1 and female is 0
                survAS=runcph(coxdf, dur, eve)
else:
    print("No Survival!")
    

