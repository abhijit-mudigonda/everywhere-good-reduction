
RIF = RealIntervalField(53)

def m_d(d):
    assert(d.is_squarefree())
    if d % 2 == 1:
        return d.abs()
    else:
        return 4*d.abs()

def eps_d(d):
    assert(d.is_squarefree())
    if d % 4 == 1:
        return 1
    return -1
    
def chi_d(d,n):
    assert(n >= 1)
    #assert(n.is_sauarefree())
    assert(d.is_squarefree())
    
    epsd = eps_d(d)
    
    result = 1
    for q,e in n.factor():
        if q % 2 == 1:
            result *= kronecker_symbol(epsd*d,q) ^ e
        else:
            if d % 2 == 1:
                result *= 1 ^ e
            else:
                result *= 0 ^ e
    return result
    
def chi01(n):
    if n % 2 == 0:
        return 0
    else:
        values = {1: 1, 3: 1, 5: -1, 7: -1}
        return values[n % 8]

def chi11(n):
    if n % 2 == 0:
        return 0
    else:
        values = {1: 1, 3: -1, 5: -1, 7: 1}
        return values[n % 8]

def chi_d_lemma44(d,n):
    assert(d.is_squarefree())
    assert(n.is_squarefree())
    if d % 2 == 1:
        ps = d.prime_factors()
        if n % 2 == 0:
            n = ZZ(n/2)
        result = prod(kronecker_symbol(n,p) for p in ps)
        return result
    else:
        ps_odd = d.odd_part().prime_factors()
        if d % 8 == 2:
            result = chi01(n)
        else:
            result = chi11(n)
        result *= prod(kronecker_symbol(n,p) for p in ps_odd)
        return result
        
def test_lemma44():
    for d in [1..1000] + [-1000..-1]:
        if d == 0:
            continue
        if not d.is_squarefree():
            continue
        for n in [1..1000]:
            if n == 0:
                continue
            if not n.is_squarefree():
                continue
            print("d,n,chi_d,chi_d_lemma44:",d,n,chi_d(d,n),chi_d_lemma44(d,n))
            assert(chi_d(d,n) == chi_d_lemma44(d,n))

#test_lemma44()

def L_chi_d(d,s):
    assert(d.is_squarefree())
    result = 1
    for p in prime_range(100):
        result *= (1-chi_d_lemma44(d,p)*p^(-s))^(-1)
    return result
    
def c_d_prime(d):
    result = 1
    result /= gamma(RIF(1/2))
    result /= RIF(zeta(2))^(1/2)
    result *= L_chi_d(d,s=RIF(1))^(1/2)
    for p in d.prime_divisors():
        result *= RIF(1+1/p)^(-1/2)
    for p in prime_range(100):
        if d % p != 0:
            if chi_d_lemma44(d,p) == 1:
                result *= RIF(1-1/p^2)^(1/2)
    return result
