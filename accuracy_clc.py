import numpy as np
import argparse
from sklearn.metrics import accuracy_score
import pandas as pd

parser = argparse.ArgumentParser(description='AMPS')

parser.add_argument('-pr', '--y_predicted', help='address to the predicted binary vector file', required=True)
parser.add_argument('-te', '--y_true', help='address to true methylation status binary vector file', required=False)
parser.add_argument('-m', '--methylation_file', help='address to true methylation file address', required=False)

args = parser.parse_args()
if args.test == None and args.methylation_file == None:
    print('Error; Provide test output')
    exit()
te = np.loadtxt(args.test)

if args.test == None:
    methylations = pd.read_table(args.methylation_file, header=None)
    methylations.columns = ['chr', 'position', 'strand', 'meth', 'unmeth', 'context', 'three']
    pr = np.asarray((methylations['meth'] / (methylations['meth'] + methylations['unmeth'])).fillna(0)).round()
else:
    pr = np.loadtxt(args.predicted)

print(accuracy_score(te.round(), pr.round()))
