import csv
import pandas as pd
import numpy as np
import random
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import pearsonr
from lifelines import CoxPHFitter




outdir = 'C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/Replication Scripts/repOut'

def runPCA(kints, df2):
    scaler = StandardScaler()
    exbic = df2[kints]
    scaler.fit(exbic)
    pcbic = scaler.transform(exbic)
    pca = PCA(n_components=1)
    pca.fit(pcbic)
    pcbic = pca.transform(pcbic) # pcbic is the principal components of the bic for each patient
    pcbic = pd.Series([i for j in range(0,len(pcbic)) for i in pcbic[j]])
    pcbic.index = df2.index
    # Make sure aligned with expression data
    muExbic = exbic.median(axis=1)
    c1 = pearsonr(muExbic, pcbic)
    if float(c1[0])<float(0):
        pcbic = -pcbic
    return pca.explained_variance_ratio_[0], pcbic


####################################################################################################################
###
### Replicate Biclusters 
###
####################################################################################################################

# Read in genes for each cluster
df1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/Replication Scripts/postProcessed_ACC_pita.csv')
# Bicluster is the index and lists of genes in the bicluster is the value
biclustMembership = pd.DataFrame(df1["ENTREZ"].apply(lambda x: [int(i) for i in x.split()]))
biclustMembership = biclustMembership.set_index(df1['Bicluster'])

#Read in a second dataset (REPLICATION DATASET
ratSec = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed exprs/GSE49278_exprs_entrez.csv',header=0, index_col=0)

#mask=[]
# mask is false if sum is zero and true if sum is not zero
#mask=[sum(ratSec.loc[i])!=0 for i in ratSec.index] 
#ratSec=ratSec[mask]

#InbMentrez=pd.Index([j for i in biclustMembership['ENTREZ'] for j in i])
d = []
intersection = pd.Series(d)
for i in np.arange(1,biclustMembership.shape[0]+1):
    intersection[i] = np.intersect1d(biclustMembership['ENTREZ'][i],ratSec.index.values)

p1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed pData/GSE49278_pData_recon.csv',header=0,index_col=0)

ks = biclustMembership.shape[0]
outNames = ['n.rows','overlap.rows','pc1.var.exp','avg.pc1.var.exp','pc1.perm.p','survival','survival.p','survival.age','survival.age.p','survival.age.sex','survival.age.sex.p']
permutations = 100
repOut = pd.DataFrame(columns=outNames)
nrows = pd.Series(biclustMembership.ENTREZ.apply(lambda x: len(x)))
orows = intersection.apply(lambda x: len(x))
p1ve = []
av_ve = []
p1pp = []
pverandoDF = pd.DataFrame([])
df2 = ratSec.T
iters = 100 #len(intersection)
for k in np.arange(0,iters):
    kints = intersection[k+1]
    testEmrows = pd.Series([])
    testEmrows[1] = kints
    pverando = []
    if len(kints) > 1: 
        varex, pcbic = runPCA(kints, df2)
        p1ve.append(varex)
        for i in np.arange(2,permutations+2):
            testEmrows[i] = random.sample(list(ratSec.index),orows[k+1])
            #need only exprs for the entrez in testEmrows
            varex, pcrando = runPCA(testEmrows[i],df2)
            pverando.append(varex)
        av_ve.append(sum(pverando)/len(pverando))
        p1pp.append(sum(pverando > varex)/permutations)
        #pverandoS=pd.Series(np.array(pverando))
        #pverandoS.columns= 'bic' +str(k)     
        #pverandoDF=pd.concat([pverandoDF,pverandoS], axis=1, sort=False)

        #############################################################################################################
        ###
        ### Replicate survival
        ###
        #############################################################################################################
        # Run if there is any survival info availible
        if any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)) and any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)) or any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)) and any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)):
            # Run if the info availible is Overall Survival
            if any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)) and any(pd.isnull(p1['OS.STATUS']).apply(lambda x: not x)):
                coxdf.append(p1['OS.STATUS'])
                coxdf.append(p1['OS.TIME'])
                dur = 'OS.TIME'
                eve = 'OS.STATUS' 
                survdf=runcph(coxdf, dur, eve)
            # If there is no Overall Survival, use Progression Free Survival
            else:
                coxdf.append(p1['PFS.STATUS'])
                coxdf.append(p1['PFS.TIME'])
                dur = 'PFS.TIME'
                eve = 'PFS.STATUS'
                survdf=runcph(coxdf, dur, eve)
            # Run if there is also age information
            if any(pd.isnull(p1['AGE']).apply(lambda x: not x)):
                coxdf.append(p1['AGE'])
            survAdf=runcph(coxdf, dur, eve)
            # Run if there is also sex information
            if any(pd.isnull(p1['SEX']).apply(lambda x: not x)):
                coxdf.append(pd.get_dummies(p1['SEX'])['M']) # Male is 1 and female is 0
                survASdf=runcph(coxdf, dur, eve)
            else:
                print("No Survival for GSE: !!!")
        surv=[surv]
##############################################################################################################################
# End Survival
##############################################################################################################################

    else:
        p1ve.append(np.NaN)
        av_ve.append(np.NaN)
        #pverando.append(np.NaN)
        p1pp.append(np.NaN)
#p1pp is the probability of a random variance explained (pve[2:len(pve)) being larger than pc1 (pve[0])
p1ve = pd.Series(p1ve, index=range(1,iters+1))
av_ve = pd.Series(av_ve, index=range(1,iters+1))
p1pp = pd.Series(p1pp, index=range(1,iters+1))

repOut = pd.DataFrame([nrows, orows, p1ve, av_ve, p1pp],  index=outNames[:5]) #, sur, surp, sur_age, surp_age, sur_agesex, surp_agesex,columns=outNames)


print(repOut)
#repOut.to_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts')


