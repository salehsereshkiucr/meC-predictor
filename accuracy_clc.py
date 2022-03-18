import numpy as np
import argparse
from sklearn.metrics import accuracy_score

parser = argparse.ArgumentParser(description='This is a demo script by nixCraft.')

parser.add_argument('-pr', '--predicted', help='address to the predicted binary vector file', required=True)
parser.add_argument('-te', '--test', help='address to test binary vector file', required=True)

args = parser.parse_args()

te = np.loadtxt(args.test)
pr = np.loadtxt(args.pr)

print(accuracy_score(te.round(), pr.round()))
