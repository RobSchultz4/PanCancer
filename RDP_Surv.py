####################################################################################################################
###
### Replicate survival
###
####################################################################################################################

#This will be the survival code implemented into main when it is working and generalized for now, I will try to get it to work on one case.

from lifelines import CoxPHFitter
from lifelines.statistics import pairwise_logrank_test
import csv
import pandas as pd
import numpy as np


p1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed pData/GSE49278_pData_recon.csv',header = 0,index_col = 0)
coxdf = []


# Run if there is any survival info availible
if any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)) and any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)) or any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)) and any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)):
    # Run if the info availible is Overall Survival
    if any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)) and any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)):
        coxdf.append(p1['OS.STATUS'])
        coxdf.append(p1['OS.TIME'])
        dur = 'OS.TIME'
        eve = 'OS.STATUS' 
    # If there is no Overall Survival, use Progression Free Survival
    else:
        coxdf.append(p1['PFS.STATUS'])
        coxdf.append(p1['PFS.TIME'])
        dur = 'PFS.TIME'
        eve = 'PFS.STATUS'
    # Run if there is also age information
    if any(pd.isnull(p1['AGE']).apply(lambda x: not x)):
        coxdf.append(p1['AGE'])
        # Run if there is also sex information
        if any(pd.isnull(p1['SEX']).apply(lambda x: not x)):
            coxdf.append(pd.get_dummies(p1['SEX'])['M']) # Male is 1 and female is 0
    coxdf = pd.DataFrame(coxdf).T
else:
    print("No Survival!")
    
cph = CoxPHFitter()
cph.fit(coxdf, duration_col = dur, event_col = eve, show_progress=True)
cph.print_summary()
