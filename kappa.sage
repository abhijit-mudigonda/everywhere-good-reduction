from tqdm import tqdm 
import argh

"""
    Run from command line via
    `sage kappa.sage $P`
    where $P is your value of P. 
"""

def kappa(P):
    """
        P: int

        Computes bounds on kappa_f for 
        f(x) = x^3-1728.
    """

    P = Integer(P)
    prod = 1
    for p in tqdm(primes(4,P)):
        if p % 3 == 1:
            leg = 1
        else:
            leg = -1
        prod *= (1 + ((2 + leg) * ((1 - pow(p,-2/3))/(pow(p,4/3)-1))))

    prod = N(prod)
    print(f"Product of factors from primes not dividing the discriminant up to {P} is:", prod)

    #factor_at_2 = (1 + ((1 - pow(2,-2/3)) * (pow(2,-1/3) + pow(2,-2/3) + 1 + pow(2,-16/3)/(1 - pow(2,-4/3)))))
    factor_at_2 = 3/4

    print("Product of local factor at 2 is:", N(factor_at_2))
    #factor_at_3 = (1 + (1-pow(3,-2/3)) * (pow(3,-1/3) + pow(3,-8/3)/(1-pow(3,-4/3))))
    factor_at_3 = 8/9 * (1 + 1/9 * pow((1 - pow(3,-2/3)),-1))

    print("Product of local factor at 3 is:",N(factor_at_3))
    prod *= factor_at_2 * factor_at_3

    f(x) = prime_pi(x)/(x-1)^(7/3)
    tail = exp( 3 * prime_pi(P)/(P-1)^(4/3) + 4 * numerical_integral(f(x), P, infinity)[0] )
    print(f"Product of factors from primes exceeding {P} is at most: ", N(tail))

    print("kappa is at least", N(2*prod), "and at most", N(2*prod*tail))

    
argh.dispatch_command(kappa)
