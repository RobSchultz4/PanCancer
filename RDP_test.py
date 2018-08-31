from subprocess import *
import pandas as pd


names = ['tcga_acc_test']  #'GSE49278']
postprocs = ['ACC_pita']
for postproc in postprocs:
    for name in names:
        a = 'python -i RDP_Main.py  --name '+name+' --postproc '+postproc
        RDPproc = Popen(a, stdout = PIPE) #, cwd='C:/Users/rschult4/Dropbox (ASU)/PanCancer')
        output = RDPproc.communicate()

tcga = pd.read_csv('postProcessed_vSurv/postProcessed_ACC_pita.csv', header = 0, index_col = 0)

######################
#This chunk of code is taken directly from the Main file. It gives out the tfs in a Series of lists of lists. 
######################

tfcols = ['Up.MEME Motif1 Correlated Matches_ACC', 'Up.MEME Motif2 Correlated Matches_ACC', 'Up.WEEDER Motif1 Correlated Matches_ACC', 'Up.WEEDER Motif2 Correlated Matches_ACC', 'TFBS_DB.Correlated Matches_ACC']
print('Done.')

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
tfs = pd.Series([ss[i].dropna().tolist() for i in ss.columns], index = ss.columns)
mask = tfs.str.len() != 0
tfs = tfs[mask]

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

for bic in repOut.index.values:
    splittfs[bic] = repOut['Replicated TFs'][bic] + repOut['Failed TFs'][bic]





#make sure none are empty
rep_tfs = pd.Series([splittfs[i].dropna().tolist() for i in splittfs.columns], index = splittfs.columns)
rep_mask = rep_tfs.str.len() != 0
rep_tfs = rep_tfs[rep_mask]


justtfs = pd.Series(index=rep_tfs.index.values)
for bic in rep_tfs.index.values:
  for tf in range(0, rep_tfs.str.len()[bic]):
  	justtfs[bic].append(rep_tfs[bic][tf])

#uniquify the tfs
rep_uassaytfs = []
for bic in rep_tfs.index.values: 
    rep_assaytfs = []
    rep_tfnums = []
    for tf in range(0, rep_tfs.str.len()[bic]):   
        rep_assaytfs.append(justtfs[bic][tf][0])
    rep_uassaytfs.append(list(set(assaytfs)))

rep_uassaytfs = pd.Series(rep_uassaytfs, index=justtfs.index.values)

##############################################
# End Chunk
##############################################

assays = ['ACC Var. Exp. First PC', 'ACC Var. Exp. First PC Perm. P-Value', 'OS_ACC','OS_ACC.p','OS.covAge_ACC','OS.covAge_ACC.p','OS.covAgeSex_ACC','OS.covAgeSex_ACC.p']
repAssay = ['pc1_var_exp', 'pc1_perm_p','os_survival', 'os_survival_p', 'os_survival_age', 'os_survival_age_p', 'os_survival_age_sex', 'os_survival_age_sex_p']
tfassays = ['tf_Pval','tf_Coef']
tests = pd.DataFrame(index = tcga['Bicluster'], columns = repAssays + tfassays)

for assay in assays:
    for bic in tcga['Bicluster']:
        tests = repOut[repAssay][bic] == tcga[assay][bic] 

for bic in tcga['Bicluster']:
    for tf in uassaytfs[bic][0]:
        tests['tf_Pval'] = uassaytfs[bic][1] == rep_uassaytfs[bic][1]
        tests['tf_Coef'] = uassaytfs[bic][1] == rep_uassaytfs[bic][1]


'''
repOut = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/PanCancerRepl/output/repOut/repOut' + name '.csv', header = 0, index_col = 0)
repstats = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/PanCancerRepl/output/repstats/repstats' + name + '.csv')
'''
# Compare the output with the tcga








