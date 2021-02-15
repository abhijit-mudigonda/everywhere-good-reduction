import timeit

@parallel
def f(n):
    return n**2

def g(n):
    return n**2

def iterator(n):
    for i in range(n):
        yield i

N = 10000

starttime = timeit.default_timer()
B = [g(n) for n in range(N)]
endtime = timeit.default_timer()
print("Time elapsed is ", endtime-starttime)

starttime = timeit.default_timer()
A = [n for n in f(iterator(N))]
endtime = timeit.default_timer()
print("Time elapsed is ", endtime-starttime)

