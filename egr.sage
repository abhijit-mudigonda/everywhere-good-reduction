#!/usr/bin/env python3

import argparse
import os
import math
from typing import Any, Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR="/home/abhijitm/research/nt/matschke_elliptic/s-unit-equations/elliptic-curve-tables/good-reduction-everywhere/"

parser = argparse.ArgumentParser()
parser.add_argument("d", type = int, action = "store", help = "the degree of extension")
args = parser.parse_args()

d = args.d

print([filename for filename in os.listdir(DATA_DIR) if filename.startswith("curves_deg_" + str(d) + "_Dmax") and filename.split('.')[1] == "sobj"])    

data_file_name = [filename for filename in os.listdir(DATA_DIR) if filename.startswith("curves_deg_" + str(d) + "_Dmax") and filename.split('.')[1] == "sobj"][-1]    
X = int(data_file_name.split('.')[0].split('_')[4][:-1])
L = load(DATA_DIR+data_file_name)

#Go through the list, marking a particular discriminant norm when you have a real (resp. imaginary) field of this 
#discriminant. At the end, generate the cumulative sum. 

discs = [[] for i in range(floor(d/2)+1)]
P.<x> = QQ[]
ctr = 0
for i in range(len(L[1])):
    ctr += 1
    K.<a> = NumberField(P(L[1][i][0][3]))
    discs[K.signature()[1]].append(abs(K.absolute_discriminant()))

sums = []
cumes = np.ndarray(shape=(floor(d/2)+1,X+1), dtype = int)
x = list(range(X+1))

for i in range(floor(d/2)+1):
    discs[i].sort()
    cume = 0
    idx = 0
    for j in range(X+1):
        while idx < len(discs[i]) and discs[i][idx] <= j:
            idx += 1
            cume += 1
        cumes[i,j] = cume
    sums.append(len(discs[i]))

    y = cumes[i,:]
    plt.plot(x,y)



#y = [0,0]+[z/(sqrt(log(z))) for z in x[2:]]
#plt.plot(x,y)

plt.savefig("./myplot.png")
print(sums)

#Want a histogram of sums and also graphs of cumulative things to see bounds


plt.bar(list(range(floor(d/2)+1)), sums, tick_label = ["("+str(d-2*s)+","+str(s)+")" for s in range(0,floor(d/2)+1)])
plt.savefig("./myhist.png")














    
    

    
    






