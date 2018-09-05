import csv
import pandas as pd
import numpy as np
import random
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import pearsonr
from lifelines import CoxPHFitter
import argparse
from multiprocessing import Pool, cpu_count, Manager

parser = argparse.ArgumentParser(description='use for testing RDP_test.py')

parser.add_argument('--name', help='Use these names to reference the material for each replication dataset', type = str)
parser.add_argument('--postproc', help='TCGA three/four letter code underscore [pita, targetscan, or tfbs_db]', type = str)
args = parser.parse_args()

def runPCA_2(kints, df):
    exbic = df[kints]
    pca = PCA(n_components=1)
    pca.fit(exbic)
    pcbic = pca.transform(exbic) # pcbic is the principal components of the bic for each patient
    pcbic = pd.Series([i for j in range(0, len(pcbic)) for i in pcbic[j]])
    pcbic.index = df.index
    # Make sure aligned with expression data
    muExbic = exbic.median(axis=1)
    c1 = pearsonr(muExbic, pcbic)
    if float(c1[0])<float(0):
        pcbic = -pcbic
    return pca.explained_variance_ratio_[0], pcbic

def repPCA_2(pcbic):
    pca = PCA(n_components=1)
    pca.fit(pcbic)
    return pca.explained_variance_ratio_[0]


#######################
# Survival Functions
#######################
def runcph(coxdf, dur, eve):
    cph = CoxPHFitter()
    cph.fit(coxdf, duration_col = dur, event_col = eve, show_progress = True, step_size = 0.5)  
    return cph

def isincluded(p1, colname):
    return any(pd.isnull(p1[colname]).apply(lambda x: not x))

def areValues(p1, colname):
    return  pd.isnull(p1[colname]).apply(lambda x: not x)

def makeMask(colname1, colname2):
    mask=[True for i in  np.arange(0, len(p1[colname1]))]
    mask1 = areValues(p1, colname1)
    mask2 = areValues(p1, colname2)
    for i in np.arange(0, len(p1[colname1])):
        if mask1[i] == False:
            mask[i] = False
        elif mask2[i] == False:
            mask[i] = False
        else:
            mask[i] = True 
    return mask

def handleNaNs(p1, colname1, colname2, coxdf):
    if any(pd.isnull(p1[colname1])) or any(pd.isnull(p1[colname2])): 
        mask = makeMask(colname1, colname2)
        coxdf = coxdf[mask]
        status = p1[colname2][mask]
        time = p1[colname1][mask]
    else:
        time = p1[colname1]
        status = p1[colname2]
    return time, status, coxdf


#################################################################################################################
###
### Setup for TF and Bicluster Correlation
###
#################################################################################################################

print('Loading postproc...')
accpp = pd.read_csv('postProcessed_vSurv/postProcessed_'+args.postproc+'.csv', header=0, index_col=0)
tfcols = ['Up.MEME Motif1 Correlated Matches_ACC', 'Up.MEME Motif2 Correlated Matches_ACC', 'Up.WEEDER Motif1 Correlated Matches_ACC', 'Up.WEEDER Motif2 Correlated Matches_ACC', 'TFBS_DB.Correlated Matches_ACC']
print('Done.')

def tfsplit(string):
    if isinstance(string, str):
        list1 = string.split(sep = ' ')
        list2 = [i.split(sep = ':') for i in list1]
    else:
        list2 = string
    return list2

inputTFs = pd.DataFrame([accpp[i] for i in tfcols]).T
# Returns a DF with rows for each bicluster and columns for each Correlated matches
# The contents are either nans for lists of lists. 
# Each sublist contains three values. The first is the entrez gene number
# The seconds is the correlation coefficient
# The third is the pvalue
#

s0 = inputTFs[inputTFs.columns[0]].apply(lambda x: tfsplit(x))
s1 = inputTFs[inputTFs.columns[1]].apply(lambda x: tfsplit(x))
s2 = inputTFs[inputTFs.columns[2]].apply(lambda x: tfsplit(x))
s3 = inputTFs[inputTFs.columns[3]].apply(lambda x: tfsplit(x))
s4 = inputTFs[inputTFs.columns[4]].apply(lambda x: tfsplit(x))

ss = pd.DataFrame([s0, s1, s2, s3, s4])

#make sure none are empty
tfs = pd.Series([ss[i].dropna().tolist() for i in ss.columns], index = ss.columns)
mask0 = tfs.str.len() != 0
tfs = tfs[mask0]

for bic in tfs.index.values:
    tfs[bic]=[i for j in tfs.loc[bic] for i in j]

#Only keep the sifnificant tfs
sigintfs = pd.Series([[] for i in range(0, len(tfs))], index=tfs.index.values)
for bic in tfs.index.values:
  for tf in range(0, tfs.str.len()[bic]):
    if float(tfs[bic][tf][2]) < 0.05:   #if there is a significant pvalue
        sigintfs[bic].append(tfs[bic][tf])

#uniquify the tfs
uassaytfs = []
for bic in tfs.index.values: 
    assaytfs = []
    tfnums = []
    for tf in range(0, tfs.str.len()[bic]):   
        assaytfs.append(sigintfs[bic][tf][0])
    uassaytfs.append(list(set(assaytfs)))

uassaytfs = pd.Series(uassaytfs, index=sigintfs.index.values)

#################################################################################################################
###
### Replicate Biclusters 
###
#################################################################################################################

# Bicluster is the index and lists of genes in the bicluster is the value
accpp.rename(columns={'Genes.1':'ENTREZ'}, inplace=True)
biclustMembership = pd.DataFrame(accpp["ENTREZ"].apply(lambda x: [int(i) for i in x.split()]))

#Read in a second dataset (REPLICATION DATASET
print('\nLoading exprs...')
ratSec = pd.read_csv('Reconstructed_exprs/' + args.name + '_exprs_entrez.csv', header=0, index_col=0)
print('Done.')


ks = biclustMembership.shape[0]
outNames = ['n_rows', 'overlap_rows', 'pc1_var_exp', 'avg_pc1_var_exp', 'pc1_perm_p', 'os_survival', 'os_survival_p', 'os_survival_age', 'os_survival_age_p', 'os_survival_age_sex', 'os_survival_age_sex_p', 'pfs_survival', 'pfs_survival_p', 'pfs_survival_age', 'pfs_survival_age_p', 'pfs_survival_age_sex', 'pfs_survival_age_sex_p', 'Replicated TFs','Failed TFs', 'Missing TFs']

df2 = ratSec.T
iters = 2 #len(intersected)
C = [] 
P = []

# Standardize the whole ratSec dataframe
scaler = StandardScaler()
exall = df2[list(set(df2.columns).intersection([i for j in biclustMembership['ENTREZ'] for i in j]))]
scaler.fit(exall)
exall_scaled= scaler.transform(exall)
df3 = pd.DataFrame(exall_scaled, columns=exall.columns, index=exall.index)

#InbMentrez=pd.Index([j for i in biclustMembership['ENTREZ'] for j in i])
intersected = pd.Series([])
for i in np.arange(1, biclustMembership.shape[0]+1):
    intersected[i] = list(set(df3.columns).intersection(biclustMembership['ENTREZ'][i]))

print('\nLoading pData...')
p1 = pd.read_csv('Reconstructed_pData/' + args.name + '_pData_recon.csv', header=0, index_col=0)
print('Done.')

import time as t1

"""def runReps(x):
    pverando.append(repPCA_2(df3[x]))
"""

permutations = 1000
repOut = pd.DataFrame(index = accpp.index.values, columns = outNames)
repOut.index.name = "Bicluster"
repOut['n_rows'] = pd.Series(accpp['Genes'])
repOut['overlap_rows'] = intersected.apply(lambda x: len(x))

for k in range(0, iters):
    start1 = t1.time()
    kints = intersected[k+1]
    pverando = []
    if len(kints) > 1: 
        start_varexp = t1.time()
        #start_pca2 = t1.time()
        varex, pcbic = runPCA_2(kints, df3)
        #end_pca2 = t1.time()
        #print('    PCA2',end_pca2-start_pca2)
        repOut.loc[k+1,'pc1_var_exp'] = varex
        #mgr = Manager()
        #pverando = mgr.list()
        rands1 = [random.sample(list(df3.columns), repOut['overlap_rows'][k+1]) for i in range(permutations)]
        #
        for x in rands1:
            #rand1 = random.sample(list(ratSec.index), repOut['overlap_rows'][k+1])
            #need only exprs for the entrez in testEmrows
            #start_rep2 = t1.time()
            varex2 = repPCA_2(df3[x])
            #end_rep2 = t1.time()
            #print('    REP2',end_rep2-start_rep2)
            pverando.append(varex2)
        #
        repOut.loc[k+1,'avg_pc1_var_exp'] = sum(pverando)/len(pverando)
        repOut.loc[k+1,'pc1_perm_p'] = sum(pverando > varex)/permutations
        print(repOut['pc1_perm_p'][k+1])
        end_varexp = t1.time()
        print('  VAREXP_time',end_varexp-start_varexp)
        #
        #####################################################################################################
        ###
        ### Replicate Survival
        ###
        #####################################################################################################
        start_surv = t1.time()
        coxdfo = pd.DataFrame([])
        coxdfp = pd.DataFrame([])
        # Run if there is any survival info availible
        if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME') or isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
            coxdfo['Bic Expression'] = pcbic
            coxdfp['Bic Expression'] = pcbic
            # Run if the info availible is Overall Survival
            if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME'):
                time, status, coxdfo = handleNaNs(p1, 'OS.TIME', 'OS.STATUS', coxdfo)
                coxdfo['OS.STATUS'] = status
                coxdfo['OS.TIME'] = time
                survo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                repOut.loc[k+1,'os_survival'] = survo.summary['z'].tolist()[0]
                repOut.loc[k+1,'os_survival_p'] = survo.summary['p'].tolist()[0]
                # Run if there is also age information
                if isincluded(p1, 'AGE'):
                    coxdfo['AGE'] = p1['AGE']
                    survAo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                    repOut.loc[k+1,'os_survival_age'] = survAo.summary['z'].tolist()[0]
                    repOut.loc[k+1,'os_survival_age_p'] = survAo.summary['p'].tolist()[0] 
                    # Run if there is also sex information
                    if isincluded(p1, 'SEX'):
                        coxdfo['SEX'] = pd.get_dummies(p1['SEX'])['M'] # Male is 1 and female is 0
                        survASo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                        repOut.loc[k+1,'os_survival_age_sex'] = survASo.summary['z'].tolist()[0]
                        repOut.loc[k+1,'os_survival_age_sex_p'] = survASo.summary['p'].tolist()[0]
            # Run if there is Progression Free Survival
            if isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
                time, status, coxdfp = handleNaNs(p1, 'PFS.TIME', 'PFS.STATUS', coxdfp)
                coxdfp['PFS.STATUS'] = status
                coxdfp['PFS.TIME'] = time 
                survp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                repOut.loc[k+1,'pfs_survival'] = survp.summary['z'].tolist()[0]   
                repOut.loc[k+1,'pfs_survival_p'] = survp.summary['p'].tolist()[0]
                # Run if there is also age information
                if isincluded(p1, 'AGE'):
                    coxdfp['AGE'] = p1['AGE']
                    survAp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                    repOut.loc[k+1,'pfs_survival_age'] = survAp.summary['z'].tolist()[0]
                    repOut.loc[k+1,'pfs_survival_age_p'] = survAp.summary['p'].tolist()[0]
                    # Run if there is also sex information
                    if isincluded(p1, 'SEX'):
                        coxdfp['SEX'] = pd.get_dummies(p1['SEX'])['M'] # Male is 1 and female is 0
                        survASp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                        repOut.loc[k+1,'pfs_survival_age_sex'] = survASp.summary['z'].tolist()[0]
                        repOut.loc[k+1,'pfs_survival_age_sex_p'] = survASp.summary['p'].tolist()[0]
        end_surv = t1.time()
        print('  SURV_time',end_surv-start_surv)
        #####################################################################################################
        ###
        ### Replicate TFs 
        ### 
        #####################################################################################################
        start_tf1 = t1.time()
        if k+1 in sigintfs.index.values: 
            reptflist = []
            failtflist = []
            tcgareptfs = list(set([int(i) for i in uassaytfs[k+1]]).intersection(df2.columns))
            repOut.loc[k+1,'Missing TFs'] = ' '.join([str(j) for j in list(set([int(i) for i in uassaytfs[k+1]]).difference(tcgareptfs))])
            repOut.loc[k+1,'Replicated TFs'] = ''
            repOut.loc[k+1,'Failed TFs'] = ''
            #
            for i in np.arange(0, len(tcgareptfs)):
                coef, pval = pearsonr(pcbic, df2[int(tcgareptfs[i])])
                C.append(coef)
                P.append(pval)
                if pval <= 0.05:
                    repOut.loc[k+1,'Replicated TFs'] = repOut['Replicated TFs'][k+1] + ' ' + str(tcgareptfs[i]) + ':' + str(coef) + ':' + str(pval)
                    #reptflist.append(tcgareptfs[i])
                else:
                    repOut.loc[k+1,'Failed TFs'] = repOut['Failed TFs'][k+1] + ' ' + str(tcgareptfs[i]) + ':' + str(coef) + ':' + str(pval)   
                    #failtflist.append(tcgareptfs[i])
            #
            repOut.loc[k+1,'Replicated TFs'] = repOut['Replicated TFs'][k+1][1:]    
            repOut.loc[k+1,'Failed TFs'] = repOut['Failed TFs'][k+1][1:]
        end_tf1 = t1.time()
        print('  TF_time',end_tf1-start_tf1)
    #
    end1 = t1.time()
    print('Bicluster #',k+1,': ',end1-start1)


#####################################################################################################
###
### Stats
###
#####################################################################################################

num0bics = sum(repOut['overlap_rows'] < 1)
num1gbics =  sum(repOut['overlap_rows'] < 2)
num2gbics =  sum(repOut['overlap_rows'] < 3)
numsmallbics = sum(repOut['overlap_rows'] < 4)
bicperRepl = sum([repOut['overlap_rows'][i] / repOut['n_rows'][i] for i in np.arange(1, len(repOut) + 1)]) / len(repOut)



if isincluded(p1, 'OS.STATUS') and isincluded(p1, 'OS.TIME'):
    perosurR = sum([repOut['os_survival_p'].dropna()[i] < .05 for i in repOut['os_survival'].dropna().index.values])/len(repOut['os_survival'].dropna())
    perosurAR = sum([repOut['os_survival_age_p'].dropna()[i] < .05 for i in repOut['os_survival_age'].dropna().index.values])/len(repOut['os_survival_age'].dropna())
    perosurASR = sum([repOut['os_survival_age_sex_p'].dropna()[i] < .05 for i in repOut['os_survival_age_sex'].dropna().index.values])/len(repOut['os_survival_age_sex'].dropna())
else:
    perosurR = 'no os survival'
    perosurAR = 'no os survival'
    perosurASR = 'no os survival'


if isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
    perpsurR = sum([repOut['pfs_survival_p'].dropna()[i] < .05 for i in repOut['pfs_survival'].dropna().index.values])/len(repOut['pfs_survival'].dropna())
    perpsurAR = sum([repOut['pfs_survival_age_p'].dropna()[i] < .05 for i in repOut['pfs_survival_age'].dropna().index.values])/len(repOut['pfs_survival_age'].dropna())
    perpsurASR = sum([repOut['pfs_survival_age_sex_p'].dropna()[i] < .05 for i in repOut['pfs_survival_age_sex'].dropna().index.values])/len(repOut['pfs_survival_age_sex'].dropna())
else:
    perpsurR = 'no pfs survival'
    perpsurAR = 'no pfs survival'
    perpsurASR = 'no pfs survival'


perReplicatedBics = sum([repOut['pc1_perm_p'].dropna()[i] < 0.05 for i in repOut['pc1_perm_p'].dropna().index.values])/len(repOut['pc1_perm_p'].dropna())


if list(set([i for j in repOut['Failed TFs'].dropna().apply(lambda x: str(x).split(sep = ' ')) for i in j]))[0] == '':
    Umissingtfs = list(set([i for j in repOut['Failed TFs'].dropna().apply(lambda x: str(x).split(sep = ' ')) for i in j]))[1:]
else:
    Umissingtfs = list(set([i for j in repOut['Failed TFs'].dropna().apply(lambda x: str(x).split(sep = ' ')) for i in j]))

numUmissingtfs = len(Umissingtfs)

av_pearsonP = sum(P)/len(P)
av_pearsonC = sum(C)/len([C])
numbicsreptf = len(list(set([i for i in repOut['Replicated TFs'].dropna().apply(lambda x: x.split(sep = ':')[0])]))[1:])

pstats = pd.Series([num0bics, num1gbics, num2gbics, numsmallbics, bicperRepl, perReplicatedBics, perosurR, perosurAR, perosurASR, perpsurR, perpsurAR, perpsurASR, numUmissingtfs, numbicsreptf] ) 

pstats.index=['number of bics with 0 genes', 'number of bics with 1 gene or less', 'number of bics with 2 genes or less', 'number of bics with 3 genes or less', 'percent of genes per bic in replication set (excludes bics with 1 or fewer genes in replication set)', 'Percent of Bics that replicate', 'Percent of bics with replicating os', 'Percent of bics with replicating os with age', 'Percent of bics with replicating survival with age and sex', 'Percent of bics replicating pfs', 'Percent of bics replicating pfs with age', 'Percent of bics replicating pfs with age and sex', 'number of tfs missing in replication dataset', 'number of bics with a replicated tf']
#####################################################################################################
###
### Output Results
###
#####################################################################################################


repOut=repOut.replace(['', np.nan], 'NA')
pstats.to_csv('output/repStats/repStats_' + args.name + '_' + args.postproc + '.csv')
repOut.to_csv('output/repOut/repOut_' + args.name + '_' + args.postproc + '.csv')

