import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='use for testing RDP_test.py')

parser.add_argument('names', help='Use these names to reference the material for each replication dataset', type = str)
parser.add_argument('o', help='Use this argument to set the output directory', type = str)
args = parser.parse_args()
name = args.names
df = pd.DataFrame([[1,2,3],[2,3,4],[3,4,5]])
df.to_csv('C:/Users/rschult4/Dropbox (ASU)/PanCancer/code/PanCancerRepl/outputtestdf.csv' + name)
print(name)




