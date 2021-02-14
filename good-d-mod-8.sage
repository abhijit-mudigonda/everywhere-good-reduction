import argh
from tqdm import tqdm
from typing import List, Set, Any, Dict
import pickle

if __name__ == "__main__":
    good_d = pickle.load(open("good-d-3000000.pickle", 'rb'))
    pos_dist = [0]*8
    neg_dist = [0]*8
    pos_tot = 0
    neg_tot = 0
    for d in tqdm(good_d):
        if d > 0:
            pos_dist[d % 8] += 1
            pos_tot += 1
        else:
            neg_dist[d % 8] += 1
            neg_tot += 1
    
    print([float(x)/pos_tot for x in pos_dist])
    print([float(x)/neg_tot for x in neg_dist])
