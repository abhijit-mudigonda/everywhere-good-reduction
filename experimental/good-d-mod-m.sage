import argh
from tqdm import tqdm
from typing import List, Set, Any, Dict
import pickle

def getDistribution(m, infile = "good-d-1000000.pickle"):
    m = int(m)
    good_d = pickle.load(open(infile, 'rb'))
    pos_dist = [0]*m
    neg_dist = [0]*m
    pos_tot = 0
    neg_tot = 0
    for d in tqdm(good_d):
        if d > 0:
            pos_dist[d % m] += 1
            pos_tot += 1
        else:
            neg_dist[d % m] += 1
            neg_tot += 1
    
    print("Positive good d:", [round(float(x)/pos_tot,3) for x in pos_dist])
    print("Negative good d:", [round(float(x)/neg_tot,3) for x in neg_dist])

argh.dispatch_commands([getDistribution])
