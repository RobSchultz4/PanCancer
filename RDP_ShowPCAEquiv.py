import csv
import pandas as pd
import numpy as np
import random
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import pearsonr

outdir = 'C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts'

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

## Load in bicluster membership   ##
## Read in genes for each cluster ##
df1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts/postProcessed_ACC_pita.csv')
# Bicluster is the index and lists of genes in the bicluster is the value
biclustMembership = pd.DataFrame(df1["ENTREZ"].apply(lambda x: [int(i) for i in x.split()]))
biclustMembership = biclustMembership.set_index(df1['Bicluster'])

# Read in an expression dataset
ratSec = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/exprs/ACC_RNAseq_medianFiltered.tsv',sep='\t',header=0, index_col=0)
#ratSec = ratSec.rename(index=str, columns={'Unnamed: 0': 'ENTREZ'})
#ratSec = ratSec.set_index('ENTREZ')

#mask = []
# mask is false if sum is zero and true if sum is not zero
#mask = [sum(ratSec.loc[i])!=0 for i in ratSec.index] 
#ratSec = ratSec[mask]

#InbMentrez=pd.Index([j for i in biclustMembership['ENTREZ'] for j in i])
d = []
intersection = pd.Series(d)
for i in np.arange(1,biclustMembership.shape[0]+1):
    intersection[i] = np.intersect1d(biclustMembership['ENTREZ'][i],ratSec.index.values)

p1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/phenotypes_panCancer.csv',header=0,index_col=0)
#p1 = p1.set_index('GSM')

eigs1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts/samplesEigs_ACC.csv',header=0,index_col=0)

for i in range(len(intersection)):   
    kints = intersection[i+1]
    df2 = ratSec.T
    ve_pc1, pc1 = runPCA(kints, df2)
    pr1 = pearsonr(pc1, eigs1.iloc[i])
    if not pr1[0]>0.999:
        print('Uh oh!', pr1, 'Cluster #', i+1, ve_pc1)

"""

k s =biclustMembership.shape[0]
outNames = ['n.rows','overlap.rows','pc1.var.exp','avg.pc1.var.exp','pc1.perm.p','survival','survival.p','survival.age','survival.age.p','survival.age.sex','survival.age.sex.p']
permutations = 10
repOut = pd.DataFrame(columns=outNames)
nrows = pd.Series(biclustMembership.ENTREZ.apply(lambda x: len(x)))
orows = intersection.apply(lambda x: len(x))
p1ve = []
pve = []
av_ve = []
p1pp = []
pveDF = pd.DataFrame([])
df2 = ratSec.T
scaler = StandardScaler()
iters = 100 #len(intersection)+1
for k in np.arange(1,iters):
    kints = intersection[k]
    if len(kints)>1:
        testEmrows = pd.Series([])
        testEmrows[1] = kints
        ve_pc1, pc1 = runPCA(kints, df2)
        pve=[]
        for i in np.arange(2,permutations+2):
            testEmrows[i]=random.sample(list(ratSec.index),orows[k])
            #need only exprs for the entrez in testEmrows
            exrando=df2[testEmrows[i]]
            scaler.fit
            pcrando=scaler.transform(exrando)
            pca=PCA(n_components=1)
            pca.fit(pcrando)
            pcrando=pca.transform(pcrando) # pcrando is the principal components of the bic for each patient
            pcrando=pd.Series([i for j in range(0,len(pcrando)) for i in pcrando[j]])
            pcrando.index=df2.index
            pve.append(pca.explained_variance_ratio_[0])
    p1ve.append(pve[0])
    av_ve.append(sum(pve[1:len(pve)])/len(pve[1:len(pve)]))
    #ppp is the probability of a random variance explained (pve[2:len(pve)) being larger than pc1 (pve[0])
    p1pp.append(sum(pve > pve[0])/permutations)
    pveS=pd.Series(np.array(pve))     
    pveS.columns='bic' + str(k)
    pveDF=pd.concat([pveDF,pveS], axis=1, sort=False)

p1ve=pd.Series(p1ve, index=range(1,iters))
av_ve=pd.Series(av_ve, index=range(1,iters))
p1pp=pd.Series(p1pp, index=range(1,iters))


repOut=pd.DataFrame([nrows, orows, p1ve, av_ve, p1pp], index=outNames[:5]) #, sur, surp, sur_age, surp_age, sur_agesex, surp_agesex,columns=outNames)


print(repOut)
#repOut.to_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts')


"""
