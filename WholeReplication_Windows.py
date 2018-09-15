from subprocess import *
import pandas as pd
import os
import csv
from multiprocessing import Pool, cpu_count

with open('dict_TCGA_Datasets.csv', mode='r') as infile:
    reader = csv.reader(infile)
    SetsnTCGA = {row[0]:(filter(None, row[1:])) for row in reader}

with open('dict_TCGA_postProc.csv', mode='r') as infile:
    reader = csv.reader(infile)
    pPnTCGA = {row[0]:list(filter(None, row[1:])) for row in reader}

def RDProc(a):
    print(' '.join(a))
    logFile = open('log/'+a[1]+'_'+a[3]+'.log','w')
    RDPproc = Popen(' '.join(a), stdout=logFile, stderr=STDOUT, shell = True)
    output = RDPproc.communicate()

def main():
    runMe = []
    for key in SetsnTCGA.keys():
        if key != 'READ/COAD':
            for name in SetsnTCGA[key]:
                for postproc in pPnTCGA[key]:
                    if not os.path.exists('output/repOut/repOut_' + name + '_' + postproc + '.csv'):
                        runMe.append(['python RDP_Main.py --name',name,'--postproc',postproc])

    cpus = 3
    # cpus = cpu_count()
    print('There are %d CPUs available.' % cpus)
    pool = Pool(processes=cpus)
    pool.map(RDProc, runMe)
    pool.close()
    pool.join()

#for a in runMe:
#    RDProc(a)

if __name__ == "__main__":
    main()

