import argh
from tqdm import tqdm
import pickle
import unittest

MAX_Q = 1000

class computeConstantForD: 
    abs_const = 1/(gamma(1/2)*zeta(2)**(1/2))

    def __init__(self, d, sgn, bound, max_q = MAX_Q):
        """
            ARGUMENTS

            d: The value of d for which we want the constant
            sgn: Whether we want real ('r') or imaginary ('i') fields
            bound: Whether we want an upper or a lower bound
            max_q: When we compute (an upper or lower bound on) 
                the product over all q s.t. chi_d(q) = 1 during the 
                computation of c'_d, we take the terms up to this value of q. 
                1000 by default, it doesn't seem to change much if it's 
                increased. 

            INSTANCE VARIABLES

            disc: The fundamental discriminant associated to this d. It will be 
                d*eps_d if d is odd (this is always 1 mod 4) or 4deps_d if 
                d is even (this is fine since d is squarefree). 
            kronecker_d: The Kronecker character of self.disc. 
            modulus: The modulus of the associated Kronecker character. 
            d_factors: A list of the factors of d. 
        """
        self.d = d
        self.sgn = sgn 
        self.bound = bound
        self.max_q = max_q
        eps_d = 1 if d%4 == 1 else -1
        self.disc = d*eps_d if d%2 == 1 else 4*d*eps_d
        self.kronecker_d = kronecker_character(self.disc)
        self.modulus = abs(self.disc)
        self.d_factors = prime_factors(d)

    def computeConstant(self):
        """
            Computes the final constant.
            This will be either c'_dc''_{R,d} or c'_dc''_{I,d}
        """
        return self.computeFirstConstant() * self.computeSecondConstant() / (abs(self.d)*2**(len(self.d_factors)))

    def computeFirstConstant(self): 
        """
            Computes c'_d for this value of d. 
        """
        return computeConstantForD.abs_const * self.computeLocalConstant() * self.computeAnnoyingProduct()

    def computeSecondConstant(self):
        """
            Computes c''_{R,d} or c''_{I,d}
        """
        if self.sgn == "r":
            if self.d % 8 == 1 or self.d % 8 == 7:
                return 11/12
            elif self.d % 8 == 3 or self.d % 8 == 5:
                if self.d > 0:
                    return 2/3
                else:
                    return 1/6
            elif self.d % 8 == 6:
                return 1/8
            elif self.d % 8 == 2 and self.d > 0:
                return 1/8
            else:
                return 0
        else:
            if self.d % 8 == 1 and self.d > 0:
                return 11/12
            elif self.d % 8 == 7 and self.d < 0:
                return 11/12
            elif self.d % 8 == 3 and self.d < 0:
                return 1/6
            elif self.d % 8 == 5 and self.d > 0:
                return 2/3
            elif (self.d % 8 == 2 or self.d % 8 == 6) and self.d < 0:
                return 1/8
            else:
                return 0

    def computeLocalConstant(self):
        """
            Compute the value of L(chi_d, 1)^{1/2}\prod_{q | m} (1+1/q)^{-1/2}
        """
        return sqrt(self.computeLocalProduct() * abs(self.evaluateLSeriesOfChiAtOne()))

    def computeLocalProduct(self):
        """
            Computes \prod_{q | m} (1+1/q)^{-1/2}
        """
        local_prod = 1.0
        for fac in self.d_factors:
            local_prod *= (fac/(fac+1))
        return local_prod

    def evaluateLSeriesOfChiAtOne(self):
        """
            Computes the value of L(1, chi_d), using a formula
            of Dirichlet for evaluating the L-series of a 
            real quadratic character at 1. We adjust if chi_d 
            is not a Kronecker symbol (i.e. when d is \pm 3 mod 8)
        """
        def evaluateLSeriesOfKroneckerAtOne():
            """
                Computes the value at 1 of the L series of the Kronecker symbol 
                with top value eps_d * d 
            """

            if self.disc < 0:
                ans = -pi*pow(self.modulus, -3/2)
                int_sum = 0
                for j in range(1,self.modulus):
                    int_sum += self.kronecker_d(j)*j
                return ans * int_sum
            else:
                ans = -pow(self.modulus, -1/2)
                int_sum = 0
                for j in range(1,self.modulus):
                    int_sum += self.kronecker_d(j)*log(sin(j*pi/self.modulus))
                return ans * int_sum
        
        if (self.d % 8 == 3 or self.d % 8 == 5):
            return (1/3)*evaluateLSeriesOfKroneckerAtOne()
        else:
            return evaluateLSeriesOfKroneckerAtOne()

    def computeAnnoyingProduct(self):
        """
            Computes a bound on the product of (1-q^{-2}) 
            over all primes q for which chi_d(q) = 1. 

            For an upper bound, it takes all the terms for
            primes up to q_max, since each term is at most 1. 

            For a lower bound, it starts with 1/zeta(2) 
            (as if all primes were in the product) and then
            removes primes up to q_max for which chi_d(q) != 1
        """

        if self.bound == 'u':
            return self.computeAnnoyingProductUB()
        else:
            return self.computeAnnoyingProductLB()

    def computeAnnoyingProductUB(self):
        prod = 1
        for p in primes(self.max_q): 
            if self.kronecker_d(p) == 1:
                prod *= (1-1/(p**2))

        if self.d % 8 == 3 or self.d % 8 == 5: 
            prod *= (3/4)

        return sqrt(prod)

    def computeAnnoyingProductLB(self):
        prod = 1/zeta(2)
        for p in primes(self.max_q): 
            if self.kronecker_d(p) != 1:
                prod *= (p**2)/(p**2-1)

        if self.d % 8 == 3 or self.d % 8 == 5: 
            prod *= (3/4)

        return sqrt(prod)

    def returnTestDict(self):
        """
            Returns a dictionary of information which the test can check
            This is probably bad software practice but whatever, seems
            cleaner than having a bunch of getters. 
        """
        return {
                "disc": self.disc,
                "modulus": self.modulus,
                "local": self.computeLocalProduct(),
                "l_series_at_one": self.evaluateLSeriesOfChiAtOne(),
                "annoying": self.computeAnnoyingProduct(),
                "second_const": self.computeSecondConstant()
                }
                    
    
def computeTail(max_d, a = 2.5, b = 0.4):
    """
        Compute an upper bound on the sum of c'_dc''_d/(|d|2^{w(d)}). 

        c'_d is at most 1/(Gamma(1/2)zeta(2)^{1/2}) * |L(1,chi_d)|^{1/2}.
        We upper bound the latter using P\'olya-Vinogradov (we could do 
        better 

        2^{w(d)} is at least 2, hence the 2 in the denominator of f(x)

    """
    f(x) = sqrt((1/2)*log(4*x)+log(log(4*x))+1/(2*sqrt(4*x)*log(4*x))+2+euler_gamma)/(2*x)
    g(x,A,B) = A*x^B*f.diff(x)
    tail = -numerical_integral(g(x,a,b),max_d,infinity)[0]
    tail *= computeConstantForD.abs_const
    print("Tail is", N(tail))
    return tail

def computeConstant(signature, bound, good_file, max_d = 10**32, max_q = MAX_Q, a = 2.5, b = 0.4):
    """
        signature: Whether we want real ('r') or imaginary ('i') fields
        bound: Whether we want an upper or a lower bound
        good_file: Where to fetch data from
        max_d: If bound is 'u', how far up to compute before we bound the tail.
        max_q: When we compute (an upper or lower bound on) 
            the product of (1-q^{-2}) over all q s.t. chi_d(q) = 1 during the 
            computation of c'_d, we take the terms up to this value of q. 
            It doesn't seem to change much if it's increased beyond 1000.
        a: The leading coefficient in our power law upper bound on number of good d up to D
        b: The exponent in our power law upper bound on number of good d up to D

        Computes an upper or lower bound on c_R or c_I, based on the values
        of signature and bound. 

    """

    assert (signature == 'r' or signature == 'i'), "Signature must be either 'r' (real) or 'i' (imaginary)"
    assert (bound == 'u' or bound == 'l'), "Bound must be either 'u' (upper) or 'l' (lower)"

    print("Finding largest d")
    with open(good_file, 'r') as goodfile:
        file_max_d = 1
        for line in goodfile:
            file_max_d = max(abs(int(line)), max_d)
    max_d = min(max_d, file_max_d)
    print("I will take d up to", max_d)

    print("Computing constant")
    total = 0
    with open(good_file, 'r') as goodfile:
        for line in goodfile:
            d = int(line)
            if abs(d) < max_d:
                total += computeConstantForD(d, signature, bound, max_q).computeConstant()
    
    if bound == 'l':
        print("Lower bound on constant is", float(total))
    else:
        total = N(total)
        print("Computed sum is", total)
        print("Upper bound on constant is", N(total+computeTail(max_d, a, b)))

class TestComputeConstantForD(unittest.TestCase): 
    """
        For now, just has a bunch of examples for which to compute the 
        various things
    """
    def test_2ru12(self):
        test_input = (2,'r','u',12)
        expected = {
            "disc": -8,
            "modulus": 8, 
            "local": 2/3,
            "l_series_at_one": -(pi/pow(8,3/2))*(1 + 3 - 5 - 7),
            "annoying": sqrt((8/9)*(120/121)),
            "second_const": 1/8,
        }
        actual = computeConstantForD(*test_input).returnTestDict()
        for key in expected.keys():
            self.assertEqual(actual[key], expected[key], key)

    def test_minus7rl14(self):
        test_input = (-7,'r','l',14)
        expected = {
            "disc": -7,
            "modulus": 7,
            "local": 7/8,
            "l_series_at_one": -(pi/pow(7,3/2))*(1 + 2 - 3 + 4 - 5 - 6),
            "annoying": sqrt((1/zeta(2))*(9/8)*(25/24)*(49/48)*(169/168)),
            "second_const": 11/12
        }

        actual = computeConstantForD(*test_input).returnTestDict()
        for key in expected.keys():
            self.assertEqual(actual[key], expected[key], key)


    def test_7iu5(self):
        test_input = (7,'i','u',5)
        expected = {
            "disc": -7,
            "modulus": 7,
            "local": 7/8,
            "l_series_at_one": -(pi/pow(7,3/2))*(1 + 2 - 3 + 4 - 5 - 6),
            "annoying": sqrt(3/4),
            "second_const": 0
        }

        actual = computeConstantForD(*test_input).returnTestDict()
        for key in expected.keys():
            self.assertEqual(actual[key], expected[key], key)


    def test_5ru10(self):
        test_input = (5,'r','u',10)
        expected = {
            "disc": 5,
            "modulus": 5,
            "local": 5/6,
            "l_series_at_one": -(1/3)*(1/pow(5,1/2))*(log(sin(pi*1/5)) - log(sin(pi*2/5)) - log(sin(pi*3/5)) + log(sin(pi*4/5))),
            "annoying": sqrt(3/4),
            "second_const": 2/3
        }
        
        actual = computeConstantForD(*test_input).returnTestDict()
        for key in expected.keys():
            self.assertEqual(actual[key], expected[key], key)


    def test_15il15(self):
        test_input = (15,'i','l',15)
        expected = {
            "disc": -15,
            "modulus": 15,
            "local": 15/24,
            "l_series_at_one": -(pi/pow(15,3/2))*(1 + 2 + 4 - 7 + 8 - 11 - 13 - 14),
            "annoying": sqrt((1/zeta(2))*(9/8)*(25/24)*(49/48)*(121/120)*(169/168)),
            "second_const": 0
        }

        actual = computeConstantForD(*test_input).returnTestDict()
        for key in expected.keys():
            self.assertEqual(actual[key], expected[key], key)

    def test_minus6ru3(self):
        test_input = (-6,'r','u',3)
        expected = {
            "disc": 24,
            "modulus": 24,
            "local": 1/2,
            "l_series_at_one": -(1/pow(24,1/2))*(log(sin(pi/24)) + log(sin(5*pi/24)) - log(sin(7*pi/24)) - log(sin(11*pi/24)) - log(sin(13*pi/24)) - log(sin(17*pi/24)) + log(sin(19*pi/24)) + log(sin(23*pi/24))),
            "annoying": 1,
            "second_const": 0
        }

        actual = computeConstantForD(*test_input).returnTestDict()
        for key in expected.keys():
            self.assertEqual(actual[key], expected[key], key)

#TO TEST, UNCOMMENT THIS BLOCK
"""
if __name__ == "__main__":
    unittest.main()
"""

argh.dispatch_commands([computeConstant])
