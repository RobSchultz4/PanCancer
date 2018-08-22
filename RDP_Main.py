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


######################## Survival Functions:
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

def handleNaNs(p1, colname1, colname2, coxdf):
    if any(pd.isnull(p1[colname1])) or any(pd.isnull(p1[colname2])): 
        mask = makeMask(colname1,colname2)
        coxdf = coxdf[mask]
        status = p1[colname2][mask]
        time = p1[colname1][mask]
    else:
        time = p1[colname1]
        status = p1[colname2]
    return time, status, coxdf

####################################################################################################################
###
### Replicate Biclusters 
###
####################################################################################################################

# Read in genes for each cluster
df1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/postProcessed_vSurv/postProcessed_ACC_pita.csv',index_col=0)
# Bicluster is the index and lists of genes in the bicluster is the value
df1.rename(columns={'Genes.1':'ENTREZ'}, inplace=True)
biclustMembership = pd.DataFrame(df1["ENTREZ"].apply(lambda x: [int(i) for i in x.split()]))

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
outNames = ['n.rows','overlap.rows','pc1.var.exp','avg.pc1.var.exp','pc1.perm.p','os.survival','os.survival.p','os.survival.age','os.survival.age.p','os.survival.age.sex','os.survival.age.sex.p','pfs.survival','pfs.survival.p','pfs.survival.age','pfs.survival.age.p','pfs.survival.age.sex','pfs.survival.age.sex.p']
permutations = 100
nrows = pd.Series(df1['Genes'])
orows = intersection.apply(lambda x: len(x))
p1ve = []
av_ve = []
p1pp = []
pverandoDF = pd.DataFrame([])
df2 = ratSec.T
iters = 3 #len(intersection)
osur =  []
osurp =  []
osur_age =  []
osurp_age =  []
osur_agesex =  []
osurp_agesex =  []
psur =  []
psurp =  []
psur_age =  []
psurp_age =  []
psur_agesex =  []
psurp_agesex =  []

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

        #####################################################################################################
        ###
        ### Replicate Survival
        ###
        #####################################################################################################

        coxdfo = pd.DataFrame([])
        coxdfp = pd.DataFrame([])
        # Run if there is any survival info availible
        if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME') or isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
            coxdfo['Bic Expression'] = pcbic
            coxdfp['Bic Expression'] = pcbic
            # Run if the info availible is Overall Survival
            if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME'):
                time, status, coxdfo = handleNaNs(p1,'OS.TIME','OS.STATUS', coxdfo)
                coxdfo['OS.STATUS'] = status
                coxdfo['OS.TIME'] = time
                survo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                osur.append(survo.summary['z'])
                osurp.append(survo.summary['p'].tolist()[0])
                # Run if there is also age information
                if isincluded(p1, 'AGE'):
                    coxdfo['AGE'] = p1['AGE']
                    survAo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                    osur_age.append(survAo.summary['z'].tolist()[0])  
                    osurp_age.append(survAo.summary['p'].tolist()[0]) 

                    # Run if there is also sex information
                    if isincluded(p1, 'SEX'):
                        coxdfo['SEX'] = pd.get_dummies(p1['SEX'])['M'] # Male is 1 and female is 0
                        survASo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                        osur_agesex.append(survASo.summary['z'].tolist()[0])
                        osurp_agesex.append(survASo.summary['p'].tolist()[0])
            else:
                osur.append(np.NaN)
                osurp.append(np.NaN)  
                osur_age.append(np.NaN)  
                osurp_age.append(np.NaN)  
                osur_agesex.append(np.NaN)
                osurp_agesex.append(np.NaN)
            # Run if there is Progression Free Survival
            if isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
                time, status, coxdfp = handleNaNs(p1,'PFS.TIME','PFS.STATUS', coxdfp)
                coxdfp['PFS.STATUS'] = status
                coxdfp['PFS.TIME'] = time 
                survp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                psur.append(survp.summary['z'].tolist()[0])   
                psurp.append(survp.summary['p'].tolist()[0])
                # Run if there is also age information
                if isincluded(p1, 'AGE'):
                    coxdfp['AGE'] = p1['AGE']
                    survAp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                    psur_age.append(survAp.summary['z'].tolist()[0])
                    psurp_age.append(survAp.summary['p'].tolist()[0])
                    # Run if there is also sex information
                    if isincluded(p1, 'SEX'):
                        coxdfp['SEX'] = pd.get_dummies(p1['SEX'])['M'] # Male is 1 and female is 0
                        survASp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                        psur_agesex.append(survASp.summary['z'].tolist()[0])
                        psurp_agesex.append(survASp.summary['p'].tolist()[0])
            else:
                psur.append(np.NaN)
                psurp.append(np.NaN)
                psur_age.append(np.NaN)
                psurp_age.append(np.NaN)
                psur_agesex.append(np.NaN)
                psurp_agesex.append(np.NaN)

        else:
            osur.append(np.NaN)
            osurp.append(np.NaN)  
            osur_age.append(np.NaN)  
            osurp_age.append(np.NaN)  
            osur_agesex.append(np.NaN)
            osurp_agesex.append(np.NaN)
            psur.append(np.NaN)
            psurp.append(np.NaN)
            psur_age.append(np.NaN)
            psurp_age.append(np.NaN)
            psur_agesex.append(np.NaN)
            psurp_agesex.append(np.NaN)
        #####################################################################################################
        ###
        ### End Survival
        ###
        #####################################################################################################

    else:
        p1ve.append(np.NaN)
        av_ve.append(np.NaN)
        #pverando.append(np.NaN)
        p1pp.append(np.NaN)
        osur.append(np.NaN)
        osurp.append(np.NaN)  
        osur_age.append(np.NaN)  
        osurp_age.append(np.NaN)  
        osur_agesex.append(np.NaN)
        osurp_agesex.append(np.NaN)
        psur.append(np.NaN)
        psurp.append(np.NaN)
        psur_age.append(np.NaN)
        psurp_age.append(np.NaN)
        psur_agesex.append(np.NaN)
        psurp_agesex.append(np.NaN)

#p1pp is the probability of a random variance explained (pve[2:len(pve)) being larger than pc1 (pve[0])
p1ve = pd.Series(p1ve, index=range(1,iters+1))
av_ve = pd.Series(av_ve, index=range(1,iters+1))
p1pp = pd.Series(p1pp, index=range(1,iters+1))

osur = pd.Series(osur, index=range(1,iters+1))
osurp =  pd.Series(osurp, index=range(1,iters+1))
osur_age =  pd.Series(osur_age, index=range(1,iters+1))
osurp_age =  pd.Series(osurp_age, index=range(1,iters+1))
osur_agesex =  pd.Series(osur_agesex, index=range(1,iters+1))
osurp_agesex =  pd.Series(osurp_agesex, index=range(1,iters+1))
psur = pd.Series(psur, index=range(1,iters+1))
psurp =  pd.Series(psurp, index=range(1,iters+1))
psur_age =  pd.Series(psur_age, index=range(1,iters+1))
psurp_age =  pd.Series(psurp_age, index=range(1,iters+1))
psur_agesex =  pd.Series(psur_agesex, index=range(1,iters+1))
psurp_agesex =  pd.Series(psurp_agesex, index=range(1,iters+1))






repOut = pd.concat([nrows, orows, p1ve, av_ve, p1pp, osur, osurp, osur_age, osurp_age, osur_agesex, osurp_agesex, psur, psurp, psur_age, psurp_age, psur_agesex, psurp_agesex], axis=1)
repOut.columns = outNames
repOut.index = df1.index.values



#print(repOut)
#repOut.to_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts')

num0bics = sum(orows < 1)
num1gbics =  sum(orows < 2)
num2gbics =  sum(orows < 3)
numsmallbics = sum(orows < 4)
bicperRepl = sum([orows[i] / nrows[i] for i in np.arange(1, len(repOut + 1))]) / len(repOut)
perReplicatedBics2gs = sum([p1ve.dropna().iloc[i] > p1pp.dropna().iloc[i] for i in np.arange(0, len(p1ve.dropna()))])/len(p1ve.dropna())
perReplicatedBics = sum([p1ve.dropna().iloc[i] > p1pp.dropna().iloc[i] for i in np.arange(0, len(p1ve.dropna()))])/len(repOut)

perosurR = sum([osurp.dropna().iloc[i] < .05 for i in np.arange(0, len(osurp.dropna()))])/len(osurp.dropna())
perpsurR = sum([psurp.dropna().iloc[i] < .05 for i in np.arange(0, len(psurp.dropna()))])/len(psurp.dropna())
perosurAR = sum([osurp_age.dropna().iloc[i] < .05 for i in np.arange(0, len(osurp_age.dropna()))])/len(osurp_age.dropna())
perpsurAR = sum([psurp_age.dropna().iloc[i] < .05 for i in np.arange(0, len(psurp_age.dropna()))])/len(psurp_age.dropna())
perosurASR = sum([osurp_agesex.dropna().iloc[i] < .05 for i in np.arange(0, len(osurp_agesex.dropna()))])/len(osurp_agesex.dropna())
perpsurASR = sum([psurp_agesex.dropna().iloc[i] < .05 for i in np.arange(0, len(psurp_agesex.dropna()))])/len(psurp_agesex.dropna())











pstats = pd.Series[num0bics, num1gbics, num2gbics, numsmallbics, bicperRepl,perReplicatedBics2gs, perReplicatedBics, perosurR, perosurAR, perosurASR, perpsurR, perpsurAR, perpsurASR]  

pstats.index=['number of bics with 0 genes','number of bics with 1 gene or less','number of bics with 2 genes or less','number of bics with 3 genes or less','percent of genes per bic in replication set (excludes bics with 1 or fewer genes in replication set)','Percent of Bics with 2 or more genes that replicate','Percent of Bics that replicate (includes all small bics)','Percent of bics with replicating os', 'Percent of bics with replicating os with age', 'Percent of bics with replicating survival with age and sec', 'Percent of bics replicating pfs', 'Percent of bics replicating pfs with age', 'Percent of bics replicating pfs with age and sex' ]




pstats.to_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Replication Scripts')
