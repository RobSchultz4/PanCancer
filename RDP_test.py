from subprocess import *
import pandas as pd

tcga = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/postProcessed_vSurv/postProcessed_ACC_pita.csv', header = 0, index_col = 0)

names = ['GSE49278', 'GSE19750', 'GSE32894']

for name in names:
    a = 'python -i RDP_Main.py ' + name
    RDPproc = Popen(a) #, cwd='C:/Users/rschult4/Dropbox (ASU)/PanCancer')
    output = RDPproc.communicate()

    '''
    repOut = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/PanCancerRepl/output/repOut/repOut' + name '.csv', header = 0, index_col = 0)
    repstats = pd.read_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/PanCancerRepl/output/repstats/repstats' + name + '.csv')
    '''
# Compare the output with the tcga








