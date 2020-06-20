#!/usr/bin/env python3

import argparse
import os
from typing import Any, Dict, List, Tuple

DATA_DIR="/home/abhijitm/research/nt/matschke_elliptic/s-unit-equations/elliptic-curve-tables/good-reduction-everywhere/"

parser = argparse.ArgumentParser()
parser.add_argument("X", type = int, action = "store", help = "the largest |Delta_K| we are interested in")
parser.add_argument("d", type = int, action = "store", help = "the degree of extension")

args = parser.parse_args()

X = args.X
d = args.d

data_file_name = [filename for filename in os.listdir(".") if filename.startswith("curves_deg_2_Dmax") and filename.split('.')[1] == "sobj"][-1]    
if X < int(data_file_name.split('.')[0].split('_')[4]):
    print("Requested value of X exceeds the largest X that has been computed for this degree. Setting X to be ", data_file_name.split('.')[0].split('_')[4], " instead")
    X = int(data_file_name.split('.')[0].split('_')[4])

L = load(data_file_name)

#Go through the list, marking a particular discriminant norm when you have a real (resp. imaginary) field of this 
#discriminant. At the end, generate the cumulative sum. 

real_norms = [0]*10000
imaginary_norms = [0
for i in len(L[1]):
    if 
    






