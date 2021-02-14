import argh
from tqdm import tqdm
from typing import List, Set, Any, Dict
import pickle

A_RANGE = 100000
D_RANGE = 6000
Q_RANGE = 1000


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


def computeGoodD(a_range = A_RANGE, d_range = D_RANGE, outfile = "") -> List[int]:
    good_d = set()
    print("Computing good d")
    for A in tqdm(range(-a_range, a_range)):
        if not ((A % 2 == 0) and (A % 16 != 0)  and (A % 16 != 4)):
            if not ((A % 3 == 0) and (A % 27 != 12)):
                d = squarefree_part(A**3-1728)
                if d == 0:
                    continue
                elif mode == 1 and ((d % 4 == 1 and sgn(d) < 0) or (d % 4 != 1 and sgn(d) > 0)):
                    continue
                else:
                    if abs(d) <= d_range:
                        good_d.add(d)
    good_d = list(good_d)
    good_d.sort()
    print(len(good_d))
    if outfile == "":
        outfile = f"good-d-{a_range}.pickle"
    with open(outfile, 'wb') as fout: 
        pickle.dump(good_d, fout)
    return good_d

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

argh.dispatch_commands([computeConstant, evaluateChi, computeGoodD, refABC, bridgeSum, testConj])
