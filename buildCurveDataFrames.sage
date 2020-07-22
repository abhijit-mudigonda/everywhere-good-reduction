#!/usr/bin/env python3

import argparse, argh
import os
import math
from typing import Any, Dict, List, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

DATA_DIR="/home/abhijitm/research/nt/matschke_elliptic/s-unit-equations/elliptic-curve-tables/good-reduction-everywhere/"
RESULT_DIR="/home/abhijitm/research/nt/matschke_elliptic/egr/code/results/"


def jInvariant(E: sage.schemes.elliptic_curves):
    if E.j_invariant().minpoly().degree() == 1:
        return E.j_invariant()
    else:
        return E.j_invariant().minpoly()

def makeRow(E: sage.schemes.elliptic_curves, K: sage.rings.number_field) -> List[Any]:
    """:
        E: An elliptic curve over K
        K: A number field
        
        Given a curve E/K, build data corresponding to the following columns:
        K: The defining polynomial for K
        d: The degree of the number field over Q
        r: The number of real places of K
        s: The number of complex places of K
        D_K: The discriminant of K
        A: The coefficient of X in a short Weierstrass model of E
        B: The constant coefficient in the same short Weierstrass model of E
        j: The j invariant of E
        j_deg: The degree of the minimal polynomial of j over Q
        gmm: Whether or not E has a global minimal model
        cm: Whether or not E has complex multiplication
    """

    row = [L[1][i][0][3], K.degree(), K.signature()[0], K.signature()[1], K.discriminant(), -K(L[1][i][1][j][0])/48, -K(L[1][i][1][j][1])/864, E.discriminant(), jInvariant(E), E.j_invariant().minpoly().degree(), E.has_global_minimal_model(), E.has_cm()]
    return row
    

def buildDfFromLoad(L: Tuple[Tuple[Any], List[Any], str]) -> pd.DataFrame:
    """
        L: The loaded .sobj file, which consists of L[0] (SageMath version), L[1] (data we want), and L[2] (runtime of computation in sage)

        returns:
            nf_at_sig: A list of floor(d/2)+1 numbers counting the number fields at each signature over which there is an elliptic curve with EGR
            curves_at_sig: A list of floor(d/2)+1 numbers counting the curves with EGR over number fields at each signature
            j_invariants: A list of floor(d/2)+1 dictionaries, each counting the number of times each j-invariant appears among EGR elliptic curves over fields of a given signature. 
    """
    P.<x> = QQ[]
    data = []
    for i in tqdm(range(len(L[1]))):
        K.<theta_K> = NumberField(P(L[1][i][0][3]))
        for j in range(len(L[1][i][1])):
            E = EllipticCurve([-K(L[1][i][1][j][0])/48, -K(L[1][i][1][j][1])/864])
            data.append(buildRow(E,K))
            df = pd.DataFrame(data, columns = ['K', 'd', 'r', 's', 'D_K', 'A','B', 'j', 'j_deg', 'gmm'])
    return df

def addColumn(d_min = 2, d_max = 5) -> None:
    """
        d_min: The lowest field extension degree of interest
        d_max: The highest field extension degree of interest

        A function for adding a column to the dataframe for each degree. 
        Edit the marked line as appropriate for the column you want to add

        WHENEVER YOU USE THIS FUNCTION, ALSO UPDATE makeRow()!
    """

    P.<x> = QQ[]
    for d in range(d_min, d_max+1):
        df = pd.read_csv(RESULT_DIR+"csv/curves_deg_"+str(d)+".csv")

        #CHANGE THIS LINE 
        df['cm'] = EllipticCurve([NumberField(P(df['K']), 'theta_K')(df['A']), NumberField(P(df['K']), 'theta_K')(df['B'])]).has_cm()

        df.to_csv(RESULT_DIR+"csv/curves_deg_"+str(d)+".csv")
    

def makeDfsFromScratch(d_min = 2, d_max = 5, datadir = DATA_DIR):
    """
        d_min: The lowest field extension degree of interest
        d_max: The highest field extension degree of interest
        datadir: Where the data should come from

        Loads the elliptic curve data with the most curves for each degree in range from the corresponding .sobj file. 
        Then calls buildDfFromLoad on each loaded .sobj to build a dataframe and make a .csv
    """
    
    for d in range(d_min, d_max+1):
        print("d = ", d)
        data_file_list = [filename for filename in os.listdir(datadir) if filename.startswith("curves_deg_" + str(d) + "_Dmax") and filename.split('.')[1] == "sobj"]    

        max_disc = 0
        for data_file in data_file_list:
            if int(data_file.split('_')[4][:-6]) > max_disc:
                data_file_name = data_file
                max_disc = int(data_file.split('_')[4][:-6])

        X = int(data_file_name.split('.')[0].split('_')[4][:-1])
        print("Using data file at: ", datadir+data_file_name)

        L = load(datadir+data_file_name)
        df = buildDfFromLoad(L)
        df.to_csv(RESULT_DIR+"csv/curves_deg_"+str(d)+".csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    argh.add_commands(parser, [makeDfsFromScratch, addColumn])
    argh.dispatch(parser)
