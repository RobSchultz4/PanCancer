import csv
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import ReplicationDatasetPermutation_noSurv as RDP

#Need to run ReplicationDatasetPermutation_noSurv.py first. Then run this:
sortedvars=pd.Series([i for i in RDP.repOut.loc['pc1.var.exp'].dropna()]).sort_values(ascending=False)


outdir = 'C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts'
dataset = 'C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts/postProcessed_ACC_pita.csv'

# Read in genes for each cluster
df1= pd.read_csv(dataset)
    # Bicluster is the index and lists of genes in the bicluster is the value
biclustMembership=pd.DataFrame(df1["ENTREZ"].apply(lambda x: [int(i) for i in x.split()]))
biclustMembership=biclustMembership.set_index(df1['Bicluster'])

#Read in a second dataset (REPLICATION DATASET
ratSec = pd.read_csv("C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed exprs/GSE49278_exprs_entrez.csv")
ratSec=ratSec.rename(index=str, columns={'Unnamed: 0': "ENTREZ"})
ratSec=ratSec.set_index('ENTREZ')

mask=[]
# mask is false if sum is zero and true if sum is not zero
mask=[sum(ratSec.loc[i])!=0 for i in ratSec.index] 
ratSec=ratSec[mask]

d=[]
intersection=pd.Series(d)
for i in np.arange(1,biclustMembership.shape[0]+1):
    intersection[i]=np.intersect1d(biclustMembership['ENTREZ'][i],ratSec.index.values)


p1 = pd.read_csv("C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed pdata/GSE49278_pData_recon.csv")
p1=p1.set_index('GSM')

ks=biclustMembership.shape[0]
outNames = ['n.rows','overlap.rows','pc1.var.exp','avg.pc1.var.exp','pc1.perm.p','survival','survival.p','survival.age','survival.age.p','survival.age.sex','survival.age.sex.p']
permutations=10
repOut=pd.DataFrame(columns=outNames)
nrows=pd.Series(biclustMembership.ENTREZ.apply(lambda x: len(x)))
orows=intersection.apply(lambda x: len(x))
pveDF=pd.DataFrame([])
df2=ratSec.T
scaler=StandardScaler()
exvar=[]
Plotbics=[i for j in [[i for i in sortedvars.index[:3]],[i for i in sortedvars.index[-3:]]] for i in j]
count=0
for k in Plotbics:
    plotdf=pd.DataFrame([])
    kints=intersection[k]
   # if len(kints)>1:
    testEmrows=pd.Series([])
    testEmrows[1]=kints
    exbic=df2[kints]
    exbic=scaler.fit_transform(exbic)
    pca=PCA(n_components=1)
    pca.fit(exbic)
    pcbic=pca.transform(exbic) # pcbic is the principal components of the bic for each patient
    pcbic=pd.Series([i for j in range(0,len(pcbic)) for i in pcbic[j]])
    pcbic.index=df2.index
    pcbic=pcbic.sort_values(ascending=True)        
    pve=[]
    for i in kints:
        exmember=df2[i]
        plotdf=pd.concat([plotdf, exmember],axis=1, sort=False)
    exvar.append(str(pca.explained_variance_ratio_))
    ind=plotdf.index
    trans=scaler.transform(plotdf)
    plotdf=pd.DataFrame(trans)
    plotdf.index=ind
    plotdf=plotdf.reindex(pcbic.index)
    #plt.figure()
    plotdf.plot.line(title='Gene Members of bic '+ str(k), legend=True, )
    plt.plot(pcbic, label='Bicluster var_exp: ' + str(exvar[count]), linewidth=5)
    plt.legend()
    count=count+1

    
plt.show()


'''
plt.plot(pcbic, label='Bicluster var_exp: ' + str(pca.explained_variance_ratio_), linewidth=5)
plt.legend()
plt.show()
'''

'''
plotdf.plot.line(title='Gene Members', legend=True, )
plt.plot(pcbic, label='Bicluster var_exp: ' + str(pca.explained_variance_ratio_), linewidth=5)
plt.legend()
plt.show()
'''
