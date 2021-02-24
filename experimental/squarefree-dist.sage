import argh
from tqdm import tqdm
from typing import List, Set, Any, Dict


def smol_count(f,R,eps):
    """
        f: The integer polynomial for which we compute f(1), f(2), ...
        R: The max value of r for which we will compute f(r)
        eps: Plot histogram with bins of width R^{1-eps}

        Plots the number of r <= x for which sqf(f(r)) < r^{1-eps} for x <= R.
    """
    R = int(R)
    eps = float(eps)
    P.<x> = QQ[]
    f = P(f)
    outlist = []
    plot_smols = []
    smols = 0
    for r in tqdm(range(1,R+1)):
        if squarefree_part(f(r)) < ceil(r**(1-eps)):
            smols += 1
            plot_smols.append((r, smols))
    plt = plot([])
    plt += list_plot(plot_smols, color = 'red')

    var('a','b')
    model(x) = a*(x**b)
    sol = find_fit(plot_smols, model)
    g(x) = model(a = sol[0].rhs(), b = sol[1].rhs())
    plt += plot(g(x), x, [1,R], color = 'blue')
    plt.save("smols.png")

def sqf_dist(f, R, eps = 1/2):
    """
        f: The integer polynomial for which we compute f(1), f(2), ...
        R: The max value of r for which we will compute f(r)
        eps: Plot histogram with bins of width R^{1-eps}
    """

    R = int(R)
    f = "x^2-108"
    P.<x> = QQ[]
    f = P(f)
    outlist = []
    for r in tqdm(range(1,R+1)):
        if f(r) != 0:
            outlist.append(squarefree_part(f(r))/f(r))
        else:
            outlist.append(0)
    h = histogram(outlist, bins = floor(R**(eps)), weights = [1/R]*R)
    h.save("smols.png")

argh.dispatch_commands([sqf_dist, smol_count])
