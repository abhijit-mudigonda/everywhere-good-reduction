#!/usr/bin/env python3

import argparse
import os
import math
from typing import Any, Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR="/home/abhijitm/research/nt/matschke_elliptic/s-unit-equations/elliptic-curve-tables/good-reduction-everywhere/"

parser = argparse.ArgumentParser()
parser.add_argument("d", type = int, action = "store", help = "The degree of the number fields I should look in")
parser.add_argument("--datadir", type = str, action = "store", default = DATA_DIR, help = "Where I can find the input .sobj file") 
parser.add_argument("--outfile", type = str, action = "store", default = "", help = "Where I should write the j invariant lists") 
parser.add_argument("--top_j", type = int, action = "store", default = 20, help = "Given K, I will print the K j invariants that occur most often for each signature")
args = parser.parse_args()

d = args.d
datadir = args.datadir
if args.outfile == "":
    outfile = "./common_j_invariants_deg_"+str(d)+".txt"
else:
    outfile = args.outfile
top_j = args.top_j

data_file_name = [filename for filename in os.listdir(datadir) if filename.startswith("curves_deg_" + str(d) + "_Dmax") and filename.split('.')[1] == "sobj"][-1]    
X = int(data_file_name.split('.')[0].split('_')[4][:-1])
L = load(datadir+data_file_name)

nf_at_sig = [0]*(floor(d/2)+1)
curves_at_sig = [0]*(floor(d/2)+1)
j_invariants = [{} for i in range(floor(d/2)+1)]

P.<x> = QQ[]
for i in range(len(L[1])):
    K.<theta_K> = NumberField(P(L[1][i][0][3]))
    nf_at_sig[K.signature()[1]] += 1

    for j in range(len(L[1][i][1])):
        E = EllipticCurve([-27*K(L[1][i][1][j][0]), -54*K(L[1][i][1][j][1])])
        j_inv = str(E.j_invariant())
        if j_inv in j_invariants[K.signature()[1]]:
            j_invariants[K.signature()[1]][j_inv] += 1
        else:
            j_invariants[K.signature()[1]][j_inv] = 1
        curves_at_sig[K.signature()[1]] += 1

plt.bar(list(range(floor(d/2)+1)), nf_at_sig, tick_label = ["("+str(d-2*s)+","+str(s)+")" for s in range(0,floor(d/2)+1)])
plt.savefig("degree_"+str(d)+"_nf_at_sig.png")
plt.bar(list(range(floor(d/2)+1)), curves_at_sig, tick_label = ["("+str(d-2*s)+","+str(s)+")" for s in range(0,floor(d/2)+1)])
plt.savefig("degree_"+str(d)+"_curves_at_sig.png")

f = open(outfile, 'w')
f.close()
f = open(outfile, 'a')
for s, j_sig in enumerate(j_invariants):
    #Print the top_j j invariants for each signature which occur most frequently
    out_list = sorted([(j, j_sig[j]) for j in j_sig.keys()], key = lambda x: x[1], reverse = True)[:top_j]
    f.write("Signature: ("+str(d-2*s)+","+str(s)+")\n")
    f.write(str(out_list))
    f.write('\n\n')
    
f.close()


