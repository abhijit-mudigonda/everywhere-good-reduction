from colorama import Fore, Back, Style
import argparse, argh

def countGoodD(N):
    N = int(N)
    D = set()
    C = set()
    max_D = 0
    for A in range(1,N):
        C.add(squarefree_part(A**3))
        if (A % 2 == 0) and ((A % 16 != 0) and (A % 16 != 4)):
            pass
        elif (A % 3 == 0) and (A % 27 != 12):
            pass
        else:
            x = squarefree_part(A**3)
            D.add(x)
            max_D = max(x, max_D)
    """ 
    for x in sorted(list(C)):
        if x in D:
            print(Fore.GREEN + str(x))
        else:
            print(Fore.RED + str(x))
    """

    count = 0
    for x in range(1, max_D):
        if x == squarefree_part(x):
            if x in D:
                print(Fore.GREEN + str(x))
            else:
                print(Fore.RED + str(x))
            count += 1
    print("Density:", float(len(D))/count)


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    argh.dispatch_command(countGoodD)


