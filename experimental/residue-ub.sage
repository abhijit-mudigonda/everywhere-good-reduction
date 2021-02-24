#!/usr/bin/env python3

import os
import math
from tqdm import tqdm



N = 500000

if __name__ == "__main__":
    ubs = []
    cum_ub_sum = 0.0
    ub_avgs = []
    count = 0
    k = 1
    
    num_factors = [0]*(N+1)
    for p in Primes():
        if p > N:
            break
        else:
            for j in range(1,ceil(N/p)):
                num_factors[p*j] += 1
        
    for m in tqdm(range(2,N)):
        if is_squarefree(m) is True:
            count += 1
            w = num_factors[m]
            ub = float((w*ln(2*m)/(e*(2**w - 1)))**(2**w-1))
            cum_ub_sum += ub
            ubs.append((m,ub))
            ub_avgs.append((m, cum_ub_sum/count))

        """
        if m > k*1500: 
            print(f"k = {k}, m = {m}")
            plt = plot([])
            plt += list_plot(ubs, color = 'red')
            save(plt, f"residue-ub-{k}.png")

            plt = plot([])
            plt += list_plot(ub_avgs, color = 'red')
            show(plt)
            save(plt, f"residue-ub-avgs-{k}.png")

            k += 1
        """

    plt = plot([])
    plt += list_plot(ubs, color = 'red')
    show(plt)
    save(plt, "residue-ub.png")

    plt = plot([])
    plt += list_plot(ub_avgs, color = 'red')
    show(plt)
    save(plt, f"residue-ub-avgs.png")


