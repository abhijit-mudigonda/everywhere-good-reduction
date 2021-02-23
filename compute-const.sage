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
BOUND_ON_C_PRIME_D = (1/4+((euler_gamma+ln(4))/log(x)) + ((euler_gamma)*ln(4)/((log(x))^2)))*(sqrt(log(x))/(2*x))
abs_const = 1/(gamma(1/2)*zeta(2)**(1/2))

class computeConstForD: 
    abs_const = 1/(gamma(1/2)*zeta(2)**(1/2))

    def __init__(self, d: int, sgn: str, q_range = Q_RANGE):
        self.d = d
        self.sgn = sgn
        self.q_range = q_range
        eps_d = 1 if d%4 == 1 else -1
        self.kronecker_d = kronecker_character(eps_d*d)
        self.conductor = d if d%2 == 1 else 4*d




    def evaluateChi(self, n: int):
        if (d % 8 == 3 or d % 8 == 5) and n % 2 == 0:
            return self.kronecker_d(n/2)
        else:
            return self.kronecker_d(n)






def evaluateChi(d: int, n: int) -> float:
    """
        Evaluate the character chi_d(n)
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
    #Compute the value of L(chi_d, 1)^{1/2}\prod_{q | m} (1+1/q)^{-1/2}
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
        local_prod *= (fac/(fac+1))

    return sqrt(local_prod * l_value)

def computeAnnoyingProduct(d: int, q_range: int) -> float:
    prod = 1.0
    for p in primes_first_n(q_range): 
        if evaluateChi(d,p) == 1:
            prod *= (1-1/(p**2))
    return prod

def computeFirstConstant(d: int, q_range: int) -> float: 
    #Compute the constant associated to the first rule of Setzer
    #This is c'_d in the paper

    abs_const = 1/(gamma(1/2)*zeta(2)**(1/2))
    local_const = computeLocalConst(d)
    annoying_const = computeAnnoyingProduct(d, q_range)**(1/2)

    return abs_const * local_const * annoying_const 

def computeSecondConstant(d: int, sign: str) -> float:
    #Compute the constant coming from the rest of the rules of Setzer
    #this is c''_d in the paper
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
        elif (d % 8 == 2 or d % 8 == 6) and d < 0:
            return 1/8
        else:
            return 0
            
def computeConstantForD(d: int, sign: str, q_range: int) -> float:
    return computeFirstConstant(d, q_range) * computeSecondConstant(d,sign)/(abs(d)*2**(len(list(factor(d)))))


def squarefreePartOfMordell(r: int) -> int:
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
                            d = squarefreePartOfMordell(r)
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
        rename("good-d.pickle",f"good-d-{i*CHUNK_SIZE}.pickle")
        sys.exit(0)

def computeTail(max_d: int, good_d_bound_a = 2.5, good_d_bound_b = 0.4):
    assert (good_d_bound_b <= 1 and good_d_bound_b > 0), "The exponent of your bound cannot exceed 1 (and really shouldn't be 1 either)"
    good_d_bound_a = float(good_d_bound_a)
    good_d_bound_b = float(good_d_bound_b)
    max_d = int(max_d)
    
    f(x) = BOUND_ON_C_PRIME_D
    g(x,A,B) = A*x^B*f.diff(x)
    tail = -numerical_integral(g(x,good_d_bound_a,good_d_bound_b),max_d,infinity)[0]
    tail *= 1/(gamma(1/2)*zeta(2)**(1/2))
    print("Tail is", tail)
    return tail

def computeConstant(sign: str, bound: str, good_file, D = 1000, q_range = Q_RANGE, good_d_bound_a = 2.5, good_d_bound_b = 0.4) -> float:
    assert (sign == 'r' or sign == 'i'), "Sign must be either 'r' (real) or 'i' (imaginary)"
    assert (bound == 'u' or bound == 'l'), "Bound must be either 'u' (upper) or 'l' (lower)"

    total = 0.0
    print("Opening file", good_file)
    with open(good_file, 'rb') as goodfile:
        D = min(D, pickle.load(goodfile))
        good_d = pickle.load(goodfile)
    
    print("Computing constant")
    for d in tqdm(good_d):
        if abs(d) < D:
            total += computeConstantForD(d,sign,q_range) 
    if bound == 'l':
        print("Lower bound on constant is", float(total))
    else:
        print("Upper bound on constant is", float(total)+computeTail(D, good_d_bound_a, good_d_bound_b))

argh.dispatch_commands([computeConstant, evaluateChi, computeGoodD, squarefreePartOfMordell, computeTail])
