from subprocess import *
import pandas as pd


with open('dict_TCGA_Datasets.csv', mode='r') as infile:
    reader = csv.reader(infile)
    with open('dict_TCGA_Datasets.csv', mode='w') as outfile:
        SetsnTCGA = {row[0]:row[1:] for row in reader}

with open('dict_TCGA_postProc.csv', mode='r') as infile:
    reader = csv.reader(infile)
    pPnTCGA = {row[0]:row[1:] for row in reader}

    with open('dict_TCGA_postProc.csv', mode='w') as outfile:
        pPnTCGA = {row[0]:row[1:] for row in reader}

for key in SetsnTCGA.keys():
    for name in SetsnTCGA[key]:
        for postproc in pPnTCGA[key]:
            if not os.path.exists('output/repOut_' + name + postproc + '.csv'):
                a = 'python -i RDP_Main.py  --name '+name+' --postproc '+postproc
                RDPproc = Popen(a) #, cwd='C:/Users/rschult4/Dropbox (ASU)/PanCancer')
                output = RDPproc.communicate()



