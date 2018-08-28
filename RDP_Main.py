import csv
import pandas as pd
import numpy as np
import random
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import pearsonr
from lifelines import CoxPHFitter
import argparse

parser = argparse.ArgumentParser(description='use for testing RDP_test.py')

parser.add_argument('name', help='Use these names to reference the material for each replication dataset', type = str)
args = parser.parse_args()
name = args.name

def runPCA(kints, df2):
    scaler = StandardScaler()
    exbic = df2[kints]
    scaler.fit(exbic)
    pcbic = scaler.transform(exbic)
    pca = PCA(n_components=1)
    pca.fit(pcbic)
    pcbic = pca.transform(pcbic) # pcbic is the principal components of the bic for each patient
    pcbic = pd.Series([i for j in range(0, len(pcbic)) for i in pcbic[j]])
    pcbic.index = df2.index
    # Make sure aligned with expression data
    muExbic = exbic.median(axis=1)
    c1 = pearsonr(muExbic, pcbic)
    if float(c1[0])<float(0):
        pcbic = -pcbic
    return pca.explained_variance_ratio_[0], pcbic


#######################
# Survival Functions
#######################
def runcph(coxdf, dur, eve):
    cph = CoxPHFitter()
    cph.fit(coxdf, duration_col = dur, event_col = eve, show_progress = False)  
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

accpp = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/postProcessed_vSurv/postProcessed_ACC_pita.csv', header=0, index_col=0)
tfcols = ['Up.MEME Motif1 Correlated Matches_ACC', 'Up.MEME Motif2 Correlated Matches_ACC', 'Up.WEEDER Motif1 Correlated Matches_ACC', 'Up.WEEDER Motif2 Correlated Matches_ACC', 'TFBS_DB.Correlated Matches_ACC']

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

alls = [s0, s1, s2, s3, s4]
ss = pd.DataFrame(alls)

#make sure none are empty
tfs = pd.Series([ss[i].dropna().tolist() for i in ss.columns], index = ss.columns)
mask = tfs.str.len() != 0
tfs = tfs[mask]

for bic in tfs.index.values:
    tfs[bic]=[i for j in tfs.loc[bic] for i in j]

#Only keep the sifnificant tfs
sigintfs = pd.Series([[] for i in range(0, len(tfs))], index=tfs.index.values)
#tfs[1][0][0, 1, 2]
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

#repOut = pd.DataFrame([], columns = ['P-Value', 'Coefficient', 'Replicated TFs', 'Failed TFs', 'Missing TFs'], index = sigintfs.index.values)

#################################################################################################################
###
### Replicate Biclusters 
###
#################################################################################################################

# Bicluster is the index and lists of genes in the bicluster is the value
accpp.rename(columns={'Genes.1':'ENTREZ'}, inplace=True)
biclustMembership = pd.DataFrame(accpp["ENTREZ"].apply(lambda x: [int(i) for i in x.split()]))

#Read in a second dataset (REPLICATION DATASET
ratSec = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed exprs/' + name + '_exprs_entrez.csv', header=0, index_col=0)

#mask=[]
# mask is false if sum is zero and true if sum is not zero
#mask=[sum(ratSec.loc[i])!=0 for i in ratSec.index] 
#ratSec=ratSec[mask]

#InbMentrez=pd.Index([j for i in biclustMembership['ENTREZ'] for j in i])
intersection = pd.Series([])
for i in np.arange(1, biclustMembership.shape[0]+1):
    intersection[i] = np.intersect1d(biclustMembership['ENTREZ'][i], ratSec.index.values)

ks = biclustMembership.shape[0]
outNames = ['n_rows', 'overlap_rows', 'pc1_var_exp', 'avg_pc1_var_exp', 'pc1_perm_p', 'os_survival', 'os_survival_p', 'os_survival_age', 'os_survival_age_p', 'os_survival_age_sex', 'os_survival_age_sex_p', 'pfs_survival', 'pfs_survival_p', 'pfs_survival_age', 'pfs_survival_age_p', 'pfs_survival_age_sex', 'pfs_survival_age_sex_p', 'Replicated TFs','Failed TFs', 'Missing TFs']

permutations = 1000
repOut = pd.DataFrame(index = accpp.index.values, columns = outNames)
repOut['n_rows'] = pd.Series(accpp['Genes'])
repOut['overlap_rows'] = intersection.apply(lambda x: len(x))

df2 = ratSec.T
iters = len(intersection)
C = []
P = []

'''
repOut['pc1_var_exp'] = []
repOut['avg_pc1_var_exp'] = []
repOut['pc1_perm_p'] = []
repOut['os_survival'] =  []
repOut['os_survival_p'] =  []
repOut['os_survival_age'] =  []
repOut['os_survival_age_p'] =  []
repOut['os_survival_age_sex'] =  []
repOut['os_survival_age_sex_p'] =  []
repOut['pfs_survival'] =  []
repOut['pfs_survival_p']=  []
repOut['pfs_survival_age']=  []
repOut['pfs_survival']=  []
repOut['pfs_survival_age_sex']=  []
repOut['pfs_survival_age_sex_p']=  []
'''

p1 = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/Reconstructed pData/' + name + '_pData_recon.csv', header=0, index_col=0)

for k in np.arange(0, iters):
    kints = intersection[k+1]
    testEmrows = pd.Series([])
    testEmrows[1] = kints
    pverando = []
    if len(kints) > 1: 
        varex, pcbic = runPCA(kints, df2)
        repOut['pc1_var_exp'].loc[k+1] = varex
        for i in np.arange(2, permutations+2):
            testEmrows[i] = random.sample(list(ratSec.index), repOut['overlap_rows'][k+1])
            #need only exprs for the entrez in testEmrows
            varex, pcrando = runPCA(testEmrows[i], df2)
            pverando.append(varex)
        repOut['avg_pc1_var_exp'].loc[k+1] = sum(pverando)/len(pverando)
        repOut['pc1_perm_p'].loc[k+1] = sum(pverando > varex)/permutations
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
                time, status, coxdfo = handleNaNs(p1, 'OS.TIME', 'OS.STATUS', coxdfo)
                coxdfo['OS.STATUS'] = status
                coxdfo['OS.TIME'] = time
                survo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                repOut['os_survival'].loc[k+1] = survo.summary['z'].tolist()[0]
                repOut['os_survival_p'].loc[k+1] = survo.summary['p'].tolist()[0]
                # Run if there is also age information
                if isincluded(p1, 'AGE'):
                    coxdfo['AGE'] = p1['AGE']
                    survAo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                    repOut['os_survival_age'].loc[k+1] = survAo.summary['z'].tolist()[0]
                    repOut['os_survival_age_p'].loc[k+1] = survAo.summary['p'].tolist()[0] 

                    # Run if there is also sex information
                    if isincluded(p1, 'SEX'):
                        coxdfo['SEX'] = pd.get_dummies(p1['SEX'])['M'] # Male is 1 and female is 0
                        survASo = runcph(coxdfo, 'OS.TIME', 'OS.STATUS')
                        repOut['os_survival_age_sex'].loc[k+1] = survASo.summary['z'].tolist()[0]
                        repOut['os_survival_age_sex_p'].loc[k+1] = survASo.summary['p'].tolist()[0]
            else:
                '''
                MAKE NA
                repOut['os_survival'].loc[k+1] = np.NaN)
                repOut['os_survival_p'].loc[k+1] = np.NaN)  
                repOut['os_survival_age'].loc[k+1] = np.NaN)  
                repOut['os_survival_sex'].loc[k+1] = np.NaN)  
                repOut['os_survival_age_sex'].loc[k+1] = np.NaN)
                repOut['os_survival_age_sex_p'].loc[k+1] = np.NaN)
                '''
            # Run if there is Progression Free Survival
            if isincluded(p1, 'PFS.STATUS') and isincluded(p1, 'PFS.TIME'):
                time, status, coxdfp = handleNaNs(p1, 'PFS.TIME', 'PFS.STATUS', coxdfp)
                coxdfp['PFS.STATUS'] = status
                coxdfp['PFS.TIME'] = time 
                survp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                repOut['pfs_survival'].loc[k+1] = survp.summary['z'].tolist()[0]   
                repOut['pfs_survival_p'].loc[k+1] = survp.summary['p'].tolist()[0]
                # Run if there is also age information
                if isincluded(p1, 'AGE'):
                    coxdfp['AGE'] = p1['AGE']
                    survAp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                    repOut['pfs_survival_age'].loc[k+1] = survAp.summary['z'].tolist()[0]
                    repOut['pfs_survival_age_p'].loc[k+1] = survAp.summary['p'].tolist()[0]
                    # Run if there is also sex information
                    if isincluded(p1, 'SEX'):
                        coxdfp['SEX'] = pd.get_dummies(p1['SEX'])['M'] # Male is 1 and female is 0
                        survASp = runcph(coxdfp, 'PFS.TIME', 'PFS.STATUS')
                        repOut['pfs_survival_age_sex'].loc[k+1] = survASp.summary['z'].tolist()[0]
                        repOut['pfs_survival_age_sex_p'].loc[k+1] = survASp.summary['p'].tolist()[0]
            else:
                '''
                repOut['pfs_survival'].loc[k+1] = np.NaN)
                repOut['pfs_survival']p.append(np.NaN)
                repOut['pfs_survival']_age.append(np.NaN)
                repOut['pfs_survival']p_age.append(np.NaN)
                repOut['pfs_survival']_agesex.append(np.NaN)
                repOut['pfs_survival']p_agesex.append(np.NaN)
                '''
        else:
            '''
            NEED to make "NA"
            repOut['os_survival'].loc[k+1] = np.NaN)
            repOut['os_survival']p.append(np.NaN)  
            repOut['os_survival']_age.append(np.NaN)  
            repOut['os_survival']p_age.append(np.NaN)  
            repOut['os_survival']_agesex.append(np.NaN)
            repOut['os_survival']p_agesex.append(np.NaN)
            repOut['pfs_survival'].loc[k+1] = np.NaN)
            repOut['pfs_survival']p.append(np.NaN)
            repOut['pfs_survival']_age.append(np.NaN)
            repOut['pfs_survival']p_age.append(np.NaN)
            repOut['pfs_survival']_agesex.append(np.NaN)
            repOut['pfs_survival']p_agesex.append(np.NaN)
            NEED to make "NA"
            '''

        #####################################################################################################
        ###
        ### End Survival 
        ###
        #####################################################################################################
    else:
        ''' Make NA

        repOut['pc1_var_exp'].loc[k+1] = np.NaN)
        repOut['avg_pc1_var_exp'].loc[k+1] = np.NaN)
        #pverando.append(np.NaN)
        repOut['pc1_perm_p'].loc[k+1] = np.NaN)
        repOut['os_survival'].loc[k+1] = np.NaN)
        repOut['os_survival']p.append(np.NaN)  
        repOut['os_survival']_age.append(np.NaN)  
        repOut['os_survival']p_age.append(np.NaN)  
        repOut['os_survival']_agesex.append(np.NaN)
        repOut['os_survival']p_agesex.append(np.NaN)
        repOut['pfs_survival'].loc[k+1] = np.NaN)
        repOut['pfs_survival']p.append(np.NaN)
        repOut['pfs_survival']_age.append(np.NaN)
        repOut['pfs_survival']p_age.append(np.NaN)
        repOut['pfs_survival']_agesex.append(np.NaN)
        repOut['pfs_survival']p_agesex.append(np.NaN)
        '''
    #####################################################################################################
    ###
    ### Replicate TFs---- THIS NEEDS AN IF STATEMENT SAYING IF BIC IS NOT IN THE INDEX OF SIGINTFS THEN
    ### SKIP IT-- that idea may be useful for survival too instead of having so many nans 
    #####################################################################################################
     
    if k+1 in sigintfs.index.values: 
        reptflist = []
        failtflist = []
        tcgareptfs = np.intersect1d(uassaytfs[k+1], df2.columns)
        repOut['Missing TFs'][k+1] = ' '.join(np.setdiff1d(np.array(uassaytfs[k+1]), tcgareptfs))
        repOut['Replicated TFs'][k+1] = ''
        repOut['Failed TFs'][k+1] = ''

        for i in np.arange(0, len(tcgareptfs)):
            coef, pval = pearsonr(pcbic, df2[int(tcgareptfs[i])])
            C.append(coef)
            P.append(pval)
            if pval <= 0.05:
                repOut['Replicated TFs'][k+1] = repOut['Replicated TFs'][k+1] + ' ' + tcgareptfs[i] + ':' + str(coef) + ':' + str(pval)
                #reptflist.append(tcgareptfs[i])
            else:
                repOut['Failed TFs'][k+1] = repOut['Failed TFs'][k+1] + ' ' + tcgareptfs[i] + ':' + str(coef) + ':' + str(pval)   
                #failtflist.append(tcgareptfs[i])

        repOut['Replicated TFs'][k+1] = repOut['Replicated TFs'][k+1][1:]    
        repOut['Failed TFs'][k+1] = repOut['Failed TFs'][k+1][1:]
        '''        
        repOut['P-Value'].loc[k+1] = P 
        repOut['Coefficient'].loc[k+1] = C
        repOut['Replicated TFs'].loc[k+1] = reptflist
        repOut['Failed TFs'].loc[k+1] = failtflist
        repOut['Missing TFs'].loc[k+1] = missingtfs[0].tolist()
        '''
#repOut['pc1_perm_p'] is the probability of a random variance explained (pve[2:len(pve)) being larger than pc1 (pve[0])

'''Shouldnt be necessary
repOut['pc1_var_exp'] = pd.Series(repOut['pc1_var_exp'], index=range(1, iters+1))
repOut['avg_pc1_var_exp'] = pd.Series(repOut['avg_pc1_var_exp'], index=range(1, iters+1))
repOut['pc1_perm_p'] = pd.Series(repOut['pc1_perm_p'], index=range(1, iters+1))

repOut['os_survival'] = pd.Series(repOut['os_survival'], index=range(1, iters+1))
repOut['os_survival']p =  pd.Series(repOut['os_survival']p, index=range(1, iters+1))
repOut['os_survival']_age =  pd.Series(repOut['os_survival']_age, index=range(1, iters+1))
repOut['os_survival']p_age =  pd.Series(repOut['os_survival']p_age, index=range(1, iters+1))
repOut['os_survival']_agesex =  pd.Series(repOut['os_survival']_agesex, index=range(1, iters+1))
repOut['os_survival']p_agesex =  pd.Series(repOut['os_survival']p_agesex, index=range(1, iters+1))
repOut['pfs_survival'] = pd.Series(repOut['pfs_survival'], index=range(1, iters+1))
repOut['pfs_survival']p =  pd.Series(repOut['pfs_survival']p, index=range(1, iters+1))
repOut['pfs_survival']_age =  pd.Series(repOut['pfs_survival']_age, index=range(1, iters+1))
repOut['pfs_survival']p_age =  pd.Series(repOut['pfs_survival']p_age, index=range(1, iters+1))
repOut['pfs_survival']_agesex =  pd.Series(repOut['pfs_survival']_agesex, index=range(1, iters+1))
repOut['pfs_survival']p_agesex =  pd.Series(repOut['pfs_survival']p_agesex, index=range(1, iters+1))

repOut = pd.concat([repOut['n_rows'], repOut['overlap_rows'], repOut['pc1_var_exp'], repOut['avg_pc1_var_exp'], repOut['pc1_perm_p'], repOut['os_survival'], repOut['os_survival_p'], repOut['os_survival_age'], repOut['os_survival_age_p'], repOut['os_survival_age_sex'], repOut['os_survival_age_sex_p'], repOut['pfs_survival'], repOut['pfs_survival']p, repOut['pfs_survival']_age, repOut['pfs_survival']p_age, repOut['pfs_survival']_agesex, repOut['pfs_survival']p_agesex], axis=1)
repOut.columns = outNames
repOut.index = df1.index.values
df1.index.values
'''

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
perReplicatedBics = sum([repOut['pc1_var_exp'].dropna()[i] > repOut['pc1_perm_p'].dropna()[i] for i in repOut['pc1_var_exp'].dropna().index.values])/len(repOut['pc1_var_exp'].dropna())

perosurR = sum([repOut['os_survival_p'].dropna()[i] < .05 for i in repOut['os_survival'].dropna().index.values])/len(repOut['os_survival'].dropna())
perpsurR = sum([repOut['pfs_survival_p'].dropna()[i] < .05 for i in repOut['pfs_survival'].dropna().index.values])/len(repOut['pfs_survival'].dropna())
perosurAR = sum([repOut['os_survival_age_p'].dropna()[i] < .05 for i in repOut['os_survival_age'].dropna().index.values])/len(repOut['os_survival_age'].dropna())
perpsurAR = sum([repOut['pfs_survival_age_p'].dropna()[i] < .05 for i in repOut['pfs_survival_age'].dropna().index.values])/len(repOut['pfs_survival_age'].dropna())
perosurASR = sum([repOut['os_survival_age_sex_p'].dropna()[i] < .05 for i in repOut['os_survival_age_sex'].dropna().index.values])/len(repOut['os_survival_age_sex'].dropna())
perpsurASR = sum([repOut['pfs_survival_age_sex_p'].dropna()[i] < .05 for i in repOut['pfs_survival_age_sex'].dropna().index.values])/len(repOut['pfs_survival_age_sex'].dropna())

if list(set([i for j in repOut['Failed TFs'].dropna().apply(lambda x: str(x).split(sep = ' ')) for i in j]))[0] == '':
    Umissingtfs = list(set([i for j in repOut['Failed TFs'].dropna().apply(lambda x: str(x).split(sep = ' ')) for i in j]))[1:]
else:
    Umissingtfs = list(set([i for j in repOut['Failed TFs'].dropna().apply(lambda x: str(x).split(sep = ' ')) for i in j]))

numUmissingtfs = len(Umissingtfs)

av_pearsonP = sum(P)/len(P)
av_pearsonC = sum(C)/len([C])
numbicsreptf = len(list(set([i for i in repOut['Replicated TFs'].dropna().apply(lambda x: x.split(sep = ':')[0])]))[1:])

pstats = pd.Series([num0bics, num1gbics, num2gbics, numsmallbics, bicperRepl, perReplicatedBics, perosurR, perosurAR, perosurASR, perpsurR, perpsurAR, perpsurASR, numUmissingtfs, numbicsreptf] ) 

pstats.index=['number of bics with 0 genes', 'number of bics with 1 gene or less', 'number of bics with 2 genes or less', 'number of bics with 3 genes or less', 'percent of genes per bic in replication set (excludes bics with 1 or fewer genes in replication set)', 'Percent of Bics replicate', 'Percent of bics with replicating os', 'Percent of bics with replicating os with age', 'Percent of bics with replicating survival with age and sec', 'Percent of bics replicating pfs', 'Percent of bics replicating pfs with age', 'Percent of bics replicating pfs with age and sex', 'number of tfs missing in replication dataset', 'number of bics with a replicated tf']
#####################################################################################################
###
### Output Results
###
#####################################################################################################


repOut=repOut.replace(['', np.nan], 'NA')
pstats.to_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/PanCancerRepl/output/repstats/repstats' + name + '.csv')
repOut.to_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/PanCancerRepl/output/repOut/repOut' + name +'.csv')


