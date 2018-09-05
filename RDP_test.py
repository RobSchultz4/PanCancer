from subprocess import *
import pandas as pd
import numpy as np


names = ['tcga_acc_test']  #'GSE49278']
postprocs = ['ACC_pita']
for postproc in postprocs:
    for name in names:
        a = 'python RDP_Main.py  --name '+name+' --postproc '+postproc
        RDPproc = Popen(a) #, cwd='C:/Users/rschult4/Dropbox (ASU)/PanCancer')
        RDPproc.communicate()

tcga = pd.read_csv('postProcessed_vSurv/postProcessed_ACC_pita.csv', header = 0, index_col = 0)
repOut = pd.read_csv('output/repOut/repOut_tcga_acc_test.csv', header = 0, index_col = 0)


######################
#This chunk of code is taken directly from the Main file. It gives out the tfs in a Series of lists of lists. 
######################

tfcols = ['Up.MEME Motif1 Correlated Matches_ACC', 'Up.MEME Motif2 Correlated Matches_ACC', 'Up.WEEDER Motif1 Correlated Matches_ACC', 'Up.WEEDER Motif2 Correlated Matches_ACC', 'TFBS_DB.Correlated Matches_ACC']

def tfsplit(string):
    if isinstance(string, str):
        list1 = string.split(sep = ' ')
        list2 = [i.split(sep = ':') for i in list1]
    else:
        list2 = string
    return list2

inputTFs = pd.DataFrame([tcga[i] for i in tfcols]).T
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
tcga_tfs = pd.Series([ss[i].dropna().tolist() for i in ss.columns], index = ss.columns)
mask = tcga_tfs.str.len() != 0
tcga_tfs = tcga_tfs[mask]

for bic in tcga_tfs.index.values:
    tcga_tfs[bic]=[i for j in tcga_tfs.loc[bic] for i in j]

#Only keep the sifnificant tfs
sigintfs = pd.Series([[] for i in range(0, len(tcga_tfs))], index=tcga_tfs.index.values)
for bic in tcga_tfs.index.values:
  for tf in range(0, tcga_tfs.str.len()[bic]):
    if float(tcga_tfs[bic][tf][2]) < 0.05:   #if there is a significant pvalue
        sigintfs[bic].append(tcga_tfs[bic][tf])

#uniquify the tfs
tcga_uassaytfs = []
for bic in tcga_tfs.index.values: 
    assaytfs = []
    tfnums = []
    for tf in range(0, tcga_tfs.str.len()[bic]):   
        assaytfs.append(sigintfs[bic][tf][0])
    tcga_uassaytfs.append(list(set(assaytfs)))

tcga_uassaytfs = pd.Series(tcga_uassaytfs, index=sigintfs.index.values)

##############################################
# This Chunk is similar to the one above but it doesn't was simplified for the repOut 
##############################################

def tfsplit(string):
    if isinstance(string, str):
        list1 = string.split(sep = ' ')
        list2 = [i.split(sep = ':') for i in list1]
    else:
        list2 = string
    return list2


splitreptfs = repOut['Replicated TFs'].apply(lambda x: tfsplit(x))
splitfailtfs = repOut['Failed TFs'].apply(lambda x: tfsplit(x))
splittfs = pd.Series(index=repOut.index.values)
for bic in repOut.index.values:
    if isinstance(splitfailtfs[bic],list) and isinstance(splitreptfs[bic],list):
        splittfs[bic] = splitreptfs[bic] + splitfailtfs[bic]
    elif isinstance(splitreptfs[bic], list):
        splittfs[bic] = splitreptfs[bic]
    elif isinstance(splitfailtfs[bic], list):
        splittfs[bic] = splitfailtfs[bic]
    else:
        splittfs[bic] = np.nan
   

#make sure none are empty
rep_tfs = pd.Series([splittfs[i] for i in splittfs.index], index = splittfs.index)
rep_mask = rep_tfs.str.len() != 0
rep_tfs = rep_tfs[rep_mask]


justtfs = pd.Series(index=rep_tfs.index.values)
for bic in rep_tfs.index.values:
    tflist = []
    if isinstance(rep_tfs[bic],list):
        for tf in range(0, len(rep_tfs[bic])):
            tflist.append(rep_tfs[bic][tf][0])
    else:
        justtfs.loc[bic] = np.nan
    justtfs.loc[bic] = tflist


#uniquify the tfs
rep_uassaytfs = pd.Series(index=justtfs.index.values)
for bic in rep_tfs.index.values: 
    rep_assaytfs = []
    if isinstance(rep_tfs[bic],list):
        for tf in range(0, len(rep_tfs[bic])):   
            rep_assaytfs.append(justtfs[bic][tf])    
        rep_uassaytfs.loc[bic] = list(set(rep_assaytfs))



##############################################
# End Chunk
##############################################

tcgaAssays = ['ACC Var. Exp. First PC', 'ACC Var. Exp. First PC Perm. P-Value', 'OS_ACC','OS_ACC.p','OS.covAge_ACC','OS.covAge_ACC.p','OS.covAgeSex_ACC','OS.covAgeSex_ACC.p']

repAssays = ['pc1_var_exp', 'pc1_perm_p','os_survival', 'os_survival_p', 'os_survival_age', 'os_survival_age_p', 'os_survival_age_sex', 'os_survival_age_sex_p']

tfAssays = ['tf_Pval','tf_Coef']


'''
#build column names:
tRcols = []
counter = 0
for assay in assays:
    for range(0,3):
        tRcols[counter] = assay + 
        counter += 1
'''
# Need to consider the tfs and the rest differently because of the extended dimensionality of the tfs
tRsubcols = ['Replication', 'TCGA', 'Same??']
tRgroups1 = [i for j in [[i,i,i] for i in repAssays] for i in j]
tRsubcols1 = [i for i in range(0,len(tRgroups1)) for i in tRsubcols]

mIndex1 = [tRgroups1,tRsubcols1]
mIndex2 = list(zip(*mIndex1))
mIndex3 = pd.MultiIndex.from_tuples(mIndex2)

testReport = pd.DataFrame(index = tcga.index.values, columns = mIndex3)

for repAssay in repAssays:
    testReport.loc[:,(repAssay,'Replication')] = repOut[repAssay] 

for assay in range(0, len(tcgaAssays)):
    testReport.loc[:,(repAssays[assay], 'TCGA')] = tcga[tcgaAssays[assay]]

for assay in range(0, len(tcgaAssays)):
    for bic in tcga.index.values:
        if repOut[repAssays[assay]][bic] != np.nan:
            testReport.loc[bic,(repAssays[assay], 'Same??')] = round(repOut[repAssays[assay]][bic],2) == round(tcga[tcgaAssays[assay]][bic],2) 
        else:
            testReport.loc[:,(repAssays[assay], 'Same??')] = np.nan

# Build output for tf tests

tRg = []
tRsc = []

for bic in tcga.index.values:
    if isinstance(rep_uassaytfs[bic],list):
        for tf in range(0, len(rep_uassaytfs[bic])):
            tRg.append(bic)
            tRsc.append(rep_uassaytfs[bic][tf]) 
    else:
        tRg.append(bic)
        tRsc.append('NA')

mIndex1 = [tRg,tRsc]
mIndex2 = list(zip(*mIndex1))
rmIndex = pd.MultiIndex.from_tuples(mIndex2)

tRsubcols = ['Replication', 'TCGA', 'Same??']
tRgroups1 = [i for j in [[i,i,i] for i in tfAssays] for i in j]
tRsubcols1 = [i for i in range(0,len(tRgroups1)) for i in tRsubcols]

mIndex1 = [tRgroups1,tRsubcols1]
mIndex2 = list(zip(*mIndex1))
cmIndex= pd.MultiIndex.from_tuples(mIndex2)

testReporttf = pd.DataFrame(index = rmIndex, columns = cmIndex)  

#fill testReporttf



# fill Replication tfs
for bic in rep_tfs.index.values:
    if not (type(rep_tfs.loc[bic])==type(float()) and np.isnan(rep_tfs.loc[bic])):
        for tf in range(0,len(rep_tfs.loc[bic])):
            for tfreport in testReporttf.loc[bic,('tf_Coef','Replication')].index.values:
                if tfreport == rep_tfs.loc[bic][tf][0]:
                    testReporttf.loc[(bic,tfreport),('tf_Coef','Replication')] =  rep_tfs.loc[bic][tf][1]
        

for bic in rep_tfs.index.values:
    if not (type(rep_tfs.loc[bic])==type(float()) and np.isnan(rep_tfs.loc[bic])):
        for tf in range(0,len(rep_tfs.loc[bic])):
            for tfreport in testReporttf.loc[bic,('tf_Pval','Replication')].index.values:
                if tfreport == rep_tfs.loc[bic][tf][0]:
                 testReporttf.loc[(bic,tfreport),('tf_Pval','Replication')] =  rep_tfs.loc[bic][tf][2]


# fill tcga tfs
for bic in tcga_tfs.index.values:
    if not (type(tcga_tfs.loc[bic])==type(float()) and np.isnan(tcga_tfs.loc[bic])):
        for tf in range(0,len(tcga_tfs.loc[bic])):
            for tfreport in testReporttf.loc[bic,('tf_Coef','TCGA')].index.values:
                if tfreport == tcga_tfs.loc[bic][tf][0]:
                 testReporttf.loc[(bic,tfreport),('tf_Coef','TCGA')] =  tcga_tfs.loc[bic][tf][1]
        

for bic in tcga_tfs.index.values:
    if not (type(tcga_tfs.loc[bic])==type(float()) and np.isnan(tcga_tfs.loc[bic])):
        for tf in range(0,len(tcga_tfs.loc[bic])):
            for tfreport in testReporttf.loc[bic].index.values:
                if tfreport == tcga_tfs.loc[bic][tf][0]:
                 testReporttf.loc[(bic,tfreport),('tf_Pval','TCGA')] =  tcga_tfs.loc[bic][tf][2]

# fill same??


for tfAssay in tfAssays:
    for bic in tcga_tfs.index.values:
        for tfreport in testReporttf.loc[bic].index.values:
          testReporttf.loc[(bic,tfreport),(tfAssay,'Same??')] = round(float(testReporttf.loc[(bic,tfreport),(tfAssay,'Replication')]),2) == round(float(testReporttf.loc[(bic,tfreport),(tfAssay,'TCGA')]),2) 

testReport.to_csv('testReport.csv')
testReporttf.to_csv('testReporttf.csv')



