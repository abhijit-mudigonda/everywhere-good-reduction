import argh
from tqdm import tqdm
from typing import List, Set, Any, Dict, TextIO
from itertools import chain
import pickle
from os import rename
import sys
R_RANGE = 100000000000000
D_RANGE = 6000
Q_RANGE = 1000
CHUNK_SIZE = 500000

def makeLowerBoundGoodD(infile: str, outfile = "lb-good-d.pickle") -> List[int]:
    """
        infile: The input pickle file from which to load, built by computeGoodD
    """
    R = []
    for char in infile:
        if char.isdigit():
            R.append(char)
    R = int(''.join(R))
    N = R/(CHUNK_SIZE)

    good_d = []
    with open(infile, 'rb') as fin, open(outfile, 'wb') as fout:
        print("Total number I need to get through is", 125)
        for i in range(100):
            print("Chunk", i)
            new_good_d = pickle.load(fin)
            good_d.extend(list(new_good_d))
        print("Making it a set")
        good_d = set(good_d)
        print("The number of distinct good d is", len(good_d))
        good_d = list(good_d).sort()
        pickle.dump(R, fout)
        pickle.dump(good_d, fout)



#@parallel(12)
def computeGoodDHelper(r: int, fout: TextIO, mode: int):
    """
        Helper function that will be parallelized
    """
    if not ((r % 2 == 0) and (r % 16 != 0)  and (r % 16 != 4)):
        if not ((r % 3 == 0) and (r % 27 != 12)):
            d = squarefreePartOfMordell(r)
            if d != 0:
                return d


def makeUpperBoundGoodD(goodfile = "good-d-3000.pickle", annfile = "annoying-d-3000.pickle", outfile = "ub-good-d.pickle") -> List[int]:
    good_d = []
    with open(goodfile, 'rb') as goodin, open(annfile, 'rb') as annin, open(outfile, 'wb') as fout: 
        while True:
            try:
                good_d += pickle.load(goodin)
            except EOFError:
                good_num = len(good_d)
                print("The number of good d is", good_num)
                break
        while True:
            try:
                good_d += pickle.load(annin)
            except EOFError:
                print("The number of annoying d is", len(good_d) - good_num)
                break
        good_d.sort()
        pickle.dump(int(3000), fout)
        pickle.dump(good_d, fout)
    return good_d
def computeGoodDParallel(r_range = R_RANGE, text_outfile = "", pickle_outfile = "") -> List[int]:
    """
        Iterates through all r in [-R_RANGE, R_RANGE] and computes the 
        set of distinct squarefree parts (see above) of r^3-1728 which appear. 
        Writes them to a file.
    """
    print("Computing good d")

    with open(text_outfile, 'w') as f_txt, open(pickle_outfile, 'wb') as f_pickle:
        try:
            good_d = computeGoodDHelper(chain(range(block_size*i, block_size*(i+1)), range(-block_size*(i+1)+1, -block_size*i+1)))
            good_d = list(set(good_d)).sort()

        except KeyboardInterrupt:
            print("You interrupted! How rude :<")
            print("I got up to r = ", r)
            writeGoodDToFileAndReturnList(good_d, text_outfile, pickle_outfile)

    return writeGoodDToFileAndReturnList(good_d, text_outfile, pickle_outfile)

def bridgeSum(d_min: int, d_max: int) -> float:
    d_min = int(d_min)
    d_max = int(d_max) 
    total = 0.0
    for d in tqdm(range(d_min, d_max+1)):
        if is_squarefree(d):
            total += 1/(d*2**(len(factor(d))))
    print(total)

def refABC(K: int, C: float):
    K = int(K)
    C = float(C)
    eps = (4*sqrt(3)/sqrt(ln(K)*ln(ln(K))))
    print("k^{1+eps} for eps = ", float(eps))
    print(f"which is 10^"+str(float(ln(K**(1+eps))/ln(10))))

def testConj(D_FLR = 2500, infile = "good-d.pickle"):
    with open(infile, 'rb') as fin:
        good_d = pickle.load(fin)
    pos_mod_eight = {k:[] for k in range(0,8)}
    neg_mod_eight = {k:[] for k in range(0,8)}
    for d in good_d:    
        if d > 0:
            pos_mod_eight[d%8].append(d)
        else:
            neg_mod_eight[8-d%8].append(d)
        
    pos_sum = {k:0.0 for k in range(0,8)}
    neg_sum = {k:0.0 for k in range(0,8)}

    D_FLR = 2500
    for k in range(0,8):
        for d in pos_mod_eight[k]:
            if d > D_FLR:
                pos_sum[k] += computeFirstConstant(d,1000)/(d*2**len(factor(d)))
        for d in neg_mod_eight[k]:
            if d > D_FLR:
                neg_sum[k] += computeFirstConstant(d,1000)/(d*2**len(factor(d)))


    print([len(pos_mod_eight[k]) for k in range(0,8)])
    print([len(neg_mod_eight[k]) for k in range(0,8)])

    print([float(pos_sum[k]/pos_sum[2]) for k in range(0,8)])
    print([float(neg_sum[k]/neg_sum[2]) for k in range(0,8)])

