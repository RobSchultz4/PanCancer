####################################################################################################################
###
### Replicate TF and Bicluster Correlation
###
####################################################################################################################
#get columns from post process
#split

import pandas as pd
import numpy as np

accpp = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/postProcessed_vSurv/postProcessed_ACC_pita.csv', header=0, index_col=0)
tfcols = ['Up.MEME Motif1 Correlated Matches_ACC', 'Up.MEME Motif2 Correlated Matches_ACC','Up.WEEDER Motif1 Correlated Matches_ACC','Up.WEEDER Motif2 Correlated Matches_ACC','TFBS_DB.Correlated Matches_ACC']

def tfsplit(string):
    if isinstance(string,str):
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


tfs = pd.Series([ss[i].dropna().tolist() for i in ss.columns], index = ss.columns)
mask = tfs.str.len() != 0
tfs = tfs[mask]

for bic in tfs.index.values:
    tfs[bic]=[i for j in tfs.loc[bic] for i in j]

sigintfs = pd.Series([[] for i in range(0, len(tfs))], index=tfs.index.values)
#tfs[1][0][0,1,2]
for bic in tfs.index.values:
  for tf in range(0,tfs.str.len()[bic]):
    if float(tfs[bic][tf][2]) < 0.05:   #if there is a significant pvalue
        sigintfs[bic].append(tfs[bic][tf])


utfs = []
for bic in tfs.index.values: 
    umbo = []
    for tf in range(0,tfs.str.len()[bic]):   
        umbo.append(sigintfs[bic][tf][0])
        tfnums = tf
    utfs.append([list(set(umbo)),tfnums])

 # gives a list of the contents of each assay. Can extract all tfs from one bic for all assays using this

#Now what?
