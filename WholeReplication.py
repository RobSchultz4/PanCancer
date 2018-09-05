from subprocess import *
import pandas as pd
import os
import csv
with open('dict_TCGA_Datasets.csv', mode='r') as infile:
    reader = csv.reader(infile)
    SetsnTCGA = {row[0]:(filter(None, row[1:])) for row in reader}

with open('dict_TCGA_postProc.csv', mode='r') as infile:
    reader = csv.reader(infile)
    pPnTCGA = {row[0]:list(filter(None, row[1:])) for row in reader}

for key in SetsnTCGA.keys():
#for key in ['ACC']:
    for name in SetsnTCGA[key]:
        for postproc in pPnTCGA[key]:
            if not os.path.exists('output/repOut/repOut_' + name + '_' + postproc + '.csv'):
                a = 'python RDP_Main.py  --name '+name+' --postproc '+postproc
                RDPproc = Popen(a, close_fds = True) #, cwd='C:/Users/rschult4/Dropbox (ASU)/PanCancer')
                RDPproc.communicate()
# end


