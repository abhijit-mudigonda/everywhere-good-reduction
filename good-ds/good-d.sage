ds = [-10000..10000]

mw_filename = 'mw_-1728d3'
good_d_filename = 'good-ds'

@parallel
def mw_basis(d, assume_GRH = False, use_Heegner_points_in_rank_1 = False):
    
    if d==131:
        use_Heegner_points_in_rank_1 = True #as magma fails there, as well as sage's E.gens().
    
    print("d:",d)
    E = EllipticCurve([0,-1728*d^3])
    E27 = EllipticCurve([0,d^3])
    print("E:",E)
    print("E27:",E27)
    
    magma_code = ''
    if assume_GRH:
        magma_code += 'SetClassGroupBounds("GRH"); \n'
    magma_code += 'E := EllipticCurve([0,%s]); \n' % (-1728*d^3,)
    magma_code += 'r1, r2 := RankBounds(E);'
    magma_code += 'use_Heegner := %s;' % ('true' if use_Heegner_points_in_rank_1 else 'false',)
    magma_code += 'if use_Heegner and (r1 eq 1) and (r2 eq 1) and (%s le 131) then' % (d,)
    magma_code += '  E27 := EllipticCurve([0,%s]); \n' % (d^3,)
    magma_code += '  reg, i := ConjecturalRegulator(E);'
    magma_code += '  reg27, i27 := ConjecturalRegulator(E27);'
    magma_code += '  "Rs", reg, reg27;'
    magma_code += 'else'
    magma_code += '  mw, a, b := Generators(E); \n'
    magma_code += '  mwSat, c := Saturation(mw); \n'
    magma_code += '  "MW", mwSat, a, b, c; \n'
    magma_code += 'end if;'
    mResult = magma.execute(magma_code)
    print("mResult:",mResult)
    
    if mResult.startswith('MW'):
        mResult = mResult[3:] #strip the "MW "
    
        mMW, abc = mResult.split('\n')
        a,b,c = [bool(x) for x in abc.split(' ')]
        
        print("mMW:",mMW)
        print("a,b,c")
        
        mMWstr = mMW.strip('[] ').split(',')
        
        print("mMW:",mMW)
        
        mwSage = [E(list(QQ(x) for x in mPstr.strip(" ()").split(" : "))) for mPstr in mMWstr]

    elif mResult.startswith('Rs'):
        mResult = mResult[3:] #strip the "MW "
        reg, reg27 = [RR(reg_str) for reg_str in mResult.split(" ")]
        print("reg,reg27:",reg,reg27)
        if reg < reg27 + 1: #+1 gives a little slack...
            pE = gp(E)
            pP = gp.ellheegner(pE)
            P = E([QQ(c) for c in pP])
            mwSage, i, reg0 = E.saturation([P])
            assert(len(mwSage) == 1)
        else:
            pE27 = gp(E27)
            pP27 = gp.ellheegner(pE27)
            P27 = E27([QQ(c) for c in pP27])
            mwSage27, i27, reg0_27 = E27.saturation([P27])

            #The following has a bug, as E is not minimal:
            #isogeny_E27_to_E = E27.isogeny(kernel=None,codomain=E,degree=3)

            #Thus we rather do:
            #Emin = E.minimal_model()
            #isogeny_E27_to_Emin = E27.isogeny(kernel=None,codomain=Emin,degree=3)
            #iso_Emin_to_E = Emin.isomorphism_to(E)
            #mwSageUnsaturated = [iso_Emin_to_E(isogeny_E27_to_Emin(Q27)) for Q27 in mwSage27]
            #Also doesn't work! (Fails to find 3-isogeny)

            #Thus we hard-code it:
            
            assert(len(mwSage27) == 1)
            x27, y27 = mwSage27[0].xy()
            a = d^3
            x,y = 4*(y27^2 + 3*a)/x27^2, 8*y27*(y27^2 - 9*a)/x27^3
            P = E(x,y) 
            mwSageUnsaturated = [P]
            mwSage, i, reg0 = E.saturation(mwSageUnsaturated)
        a,b,c = True, True, True
    else:
        raise 

    print("mwSage:",mwSage)
    mwSageNonTorsion = [P.xy() for P in mwSage if not P.is_finite_order()]
    #mw_sat = E.saturation(mwSage)
    return mwSageNonTorsion, a, b, c #, proved

def compute_mw_bases():

    global ds
    global mw_filename

    try:
        mws = load(mw_filename+".sobj")
        keys = list(mws.keys())
        for k in keys:
            if mws[k] == "NO DATA":
                mws.pop(k)
    except (FileNotFoundError, EOFError):
        mws = {}
        

    #mws = {} #Restart whole computation!

    ds = [d for d in ds if (d not in mws) and d.is_squarefree()]
    i = 0

    try:
        for result in mw_basis(ds):
            print("result:",result)
            d = result[0][0][0]
            mw = result[1]
            print("d:",d)
            print('result:',result)

            mws[d] = mw
            
            i += 1
            if i % 100 == 0:
                save(mws,mw_filename)
                save(mws,mw_filename+"_safety")
                
            
    except KeyboardInterrupt:
        pass
    save(mws,mw_filename)
    save(mws,mw_filename+"_safety")

compute_mw_bases()

def compute_good_ds_from_mw_bases():

    mws = load(mw_filename+".sobj")

    load("mordellSInteger27.sage")

    try:
        good_ds = load(good_d_filename+".sobj")
    except (FileNotFoundError, EOFError):
        good_ds = {}
    i = 0

    try:
        for d, mw_abc in mws.items():
            if d in good_ds:
                continue    
            if mw_abc == "NO DATA":
                print("NO DATA for d =",d)
                continue
            
            mw, a, b, c = mw_abc
            
            #print("d,mw:",d,mw)
            
            E = EllipticCurve([0,-1728*d^3])
            mwE = [E(P) for P in mw]
            mwInf = [P for P in mwE if not P.is_finite_order()]
            
            #E1 = EllipticCurve([0,-1728*d^3])
            #isogeny = E0.isogeny_of_degree_3(kernel=None, codomain = E1, degree = 3)
            #mw1 = [isogeny(P) for P in E0]
            
            Ps = computeSIntegralPoints(E,S=[],mwBasis=mwInf,verbose=False,bothSigns=True)
            #print("len(Ps):",len(Ps))
            #dY^2 = X^3 - 1728
            #y^2 = x^3 - 1728d^3
            #y = Y*d^2
            #x = X*d
            
            is_good = False
            
            P0s = []
            for P in Ps:
                x, y = P.xy()
                X = x/d
                Y = y/d^2
                if X in ZZ and Y in ZZ and Y != 0:
                    #d is squarefree part of (X^3-1728):
                    if X % 2 == 0:
                        if X % 16 not in [0,4]:
                            continue
                    if X % 3 == 0:
                        if X % 27 != 12:
                            continue
                    is_good = True
            
            print("d is good:",d,is_good)
            good_ds[d] = is_good
            
            i += 1
            if i % 100 == 0:
                save(good_ds,good_d_filename)
                save(good_ds,good_d_filename+"_safety")
          
    except KeyboardInterrupt:
        pass

    save(good_ds,good_d_filename)
    save(good_ds,good_d_filename+"_safety")
        
    print("good_ds:",good_ds)

compute_good_ds_from_mw_bases()
