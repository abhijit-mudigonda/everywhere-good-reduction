import argh
from tqdm import tqdm

def posnegrange(a, b):
    for x in range(a,b):
        yield x
    for x in range(-b+1, -a+1):
        yield x

def computeGoodD(d_min = 1, d_max = 3000, good_outfile = "ub-def-good-d.txt", annoying_outfile = "ub-annoying_d.txt"):
    good_d = []
    annoying_d = []
    
    def writeListsToFile():
        with open(good_outfile, 'w') as goodfile, open(annoying_outfile, 'w') as annoyingfile:
            print("Writing to good file")
            for d in tqdm(good_d):
                goodfile.write(str(d)+'\n')

            print("Writing to annoying file")
            for d in tqdm(annoying_d):
                annoyingfile.write(str(d)+'\n')

    for d in posnegrange(d_min, d_max):
        if d.is_squarefree() and d != 0:
            E = EllipticCurve([0, -1728*(d**3)])
            try:
                int_pts = E.integral_points()
            except RuntimeError:
                print(d, "is annoying")
                annoying_d.append(d)
                continue

            if len(int_pts) > 1:
                flag = false
                for pt in int_pts:
                    x = pt[0]/d
                    if x.is_integer() and x != 12 and \
                        not (x % 2 == 0 and (x % 16 == 8 or x % 16 == 12)) and \
                        not (x % 3 == 0 and x % 27 != 12):       
                        print(d, "is good")
                        good_d.append(d)
                        flag = true
                        break
                if flag == false:
                    print(d, "is bad")
                        
            else:
                print(d, "is bad")

    writeListsToFile()

argh.dispatch_commands([computeGoodD])

