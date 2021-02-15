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

def evaluateChi(d: int, n: int) -> float:
    """
        Evaluate the character
    """
    if d%4 == 1:
        eps_d = 1
    else:
        eps_d = -1
    prod = 1
    for q in [fac[0] for fac in factor(n)]:
        if q != 2:
            prod *= legendre_symbol(d*eps_d, q)
    return prod

def computeLocalConst(d: int) -> float:
    #Compute the value of L(chi_d, 1)\prod_{q | m} (1+1/q)
    if d%2 == 1:
        m = abs(d)
    else:
        m = 4*abs(d)

    l_value = 0.0
    for i in range(1, m):
        if gcd(i,m) == 1:
            l_value += evaluateChi(d,i)*psi(i/m)
    l_value *= (-1/m)


    local_prod = 1.0
    m_factors = [tup[0] for tup in list(factor(m))]
    for fac in m_factors:
        local_prod *= (1+1/fac)

    return local_prod * l_value

def computeAnnoyingProduct(d: int, q_range: int) -> float:
    prod = 1.0
    for p in primes_first_n(q_range): 
        if evaluateChi(d,p) == 1:
            prod *= (1-1/p**2)
    return prod

def computeFirstConstant(d: int, q_range: int) -> float: 
    #Compute the constant associated to the first rule of Setzer
    #This is the constant from Lemma 5.2 of the paper. 

    abs_const = 1/(gamma(1/2)*zeta(2)**(1/2))
    local_const = computeLocalConst(d)**(1/2)
    annoying_const = computeAnnoyingProduct(d, q_range)**(1/2)

    return abs_const * local_const * annoying_const 

def computeRemainingConstant(d: int, sign: str) -> float:
    if sign == "r":
        if d % 8 == 1 or d % 8 == 7:
            return 11/12
        elif d % 8 == 3 or d % 8 == 5:
            if d > 0:
                return 2/3
            else:
                return 1/6
        elif d % 8 == 6:
            return 1/8
        elif d % 8 == 2 and d > 0:
            return 1/8
        else:
            return 0
    else:
        if d % 8 == 1 and d > 0:
            return 11/12
        elif d % 8 == 7 and d < 0:
            return 11/12
        elif d % 8 == 3 and d < 0:
            return 1/6
        elif d % 8 == 5 and d > 0:
            return 2/3
        elif d % 8 == 2 and d < 0:
            return 1/8
        elif d % 8 == 6 and d < 0:
            return 1/8
        else:
            return 0
            
def computeConstantForD(d: int, sign: str, q_range: int) -> float:
    return computeFirstConstant(d, q_range) * computeRemainingConstant(d,sign)/(abs(d)*2**(len(list(factor(d)))))


def squarefreePart(r: int) -> int:
    """
        Computes the squarefree part of f(r) = r^3-1728. 
        The squarefree part of an integer n is the unique squarefree d such that n = dt^2

        Since f(r) = (r-12)(r^2+12r+144) we factor these two separately to improve runtime
    """
    F_1 = factor(r-12)
    F_2 = factor(r**2+12*r+144)
    d = 1
    e_2 = 0
    e_3 = 0
    for (pp, e) in chain(F_1,F_2):
        if pp == 2:
            e_2 += e
        elif pp == 3:
            e_3 += e
        elif e%2:
            d *= pp
    if e_2 % 2:
        d *= 2
    if e_3 % 2:
        d *= 3
    if r < 12:
        d = -d
    return d

#@parallel(12)
def computeGoodDHelper(r: int, fout: TextIO, mode: int):
    """
        Helper function that will be parallelized
    """
    if not ((r % 2 == 0) and (r % 16 != 0)  and (r % 16 != 4)):
        if not ((r % 3 == 0) and (r % 27 != 12)):
            d = squarefreePart(r)
            if d != 0:
                return d


def writeGoodDToFileAndReturnList(good_d: Set[int], text_outfile: str, pickle_outfile: str) -> List[int]:
    good_d = list(good_d)
    good_d.sort()
    print("The number of good d is: ", len(good_d))

    if text_outfile == "":
        text_outfile = f"good-d-{r_range}.txt"
    if pickle_outfile == "":
        pickle_outfile = f"good-d-{r_range}.pickle"

        print("Writing text output file")
        for d in tqdm(good_d):
            f.write(f"{d}\n")
        f.flush()
    with open(pickle_outfile, 'wb') as f:
        print("Writing pickle output file")
        pickle.dump(good_d, f)

    return good_d

def iterator(x):
    for i in range(x):
        yield i

def computeGoodD(r_range = R_RANGE):
    """
        Iterates through all r in [-r_RANGE, r_RANGE] and computes the 
        set of distinct squarefree parts (see above) of r^3-1728 which appear. 
        Writes them to a file.
    """

    print("Computing good d")
    d_count = 0 #Upper bound on the number of good d
    try:
        with open("good-d.pickle", 'ab') as f_pickle:
            for i in range(ceil(R_RANGE/CHUNK_SIZE)):
                print("i is: ", i)
                good_d = set()
                for r in tqdm(chain(range(CHUNK_SIZE*i, CHUNK_SIZE*(i+1)), range(-CHUNK_SIZE*(i+1)+1, -CHUNK_SIZE*i+1))):
                    if not ((r % 2 == 0) and (r % 16 != 0)  and (r % 16 != 4)):
                        if not ((r % 3 == 0) and (r % 27 != 12)):
                            d = squarefreePart(r)
                            if d != 0:
                                good_d.add(d)
                d_count += len(good_d)
                print("Writing to file, do not interrupt")
                pickle.dump(good_d, f_pickle)
    except KeyboardInterrupt:
        print("You interrupted! How rude :<")
        print("I got up to |r| at least", i*CHUNK_SIZE)
        print("The number of good d so far is at most", d_count+len(good_d))

        print("Dumping")
        with open("good-d.pickle", 'ab') as f_pickle:
            pickle.dump(good_d, f_pickle)
        print("Renaming")
        rename("good-d.pickle",f"nu-good-d-{i*CHUNK_SIZE}.pickle")
        sys.exit(0)


def loadGoodD(infile: str) -> List[int]:
    """
        infile: The input pickle file from which to load, built by computeGoodD
    """
    R = []
    for char in infile:
        if char.isdigit():
            R.append(char)
    R = int(''.join(R))
    N = R/(CHUNK_SIZE)

    good_d = set()
    d_count = 0
    with open(infile, 'rb') as fin:
        for i in range(N+1):
            print(i)
            new_good_d = pickle.load(fin)
            d_count += len(new_good_d)
            good_d = union(good_d, new_good_d)
    print("The number of non-distinct items I loaded is", d_count)
    print("The number of distinct good d is", len(good_d))
    good_d = list(good_d).sort()
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

def computeConstant(sign: str, infile = "good-d.pickle", q_range = Q_RANGE) -> float:
    assert(sign == 'r' or sign == 'i')
    total = 0.0
    print("Computing constant")
    with open(infile, 'rb') as fin:
        good_d = pickle.load(fin)
    for d in tqdm(good_d):
        total += computeConstantForD(d,sign,q_range) 
    print(float(total))


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

argh.dispatch_commands([computeConstant, evaluateChi, computeGoodD, refABC, bridgeSum, testConj, squarefreePart, loadGoodD])
