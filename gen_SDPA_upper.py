#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:25:15 2020

@author: pemeriau


In this file, we create a SDP and output it as a text file with
extension dat-s to it is readable by SDPA-GMP.
"""

# Import useful libraries. 
import picos as pc
from math import factorial, sqrt
import numpy as np
from scipy.special import binom
import os


def get_Qi(Q,i,const_ij,m):
    """
    Aim: 
    ----    
     Equalising two polynomials where one is obtained by a SOS 
     decomposition in the canonical basis and the other one is expressed
     in the Laguerre basis.
     
    Parameters
    ----------
     Q : matrix for the SOS decomposition
     i : integer
         degree at which we compte the coefficients.     
     const_ij : list
                contains indices of Q at which coefficients i+j= const.

    Returns
    -------
     Real that is a sum of coefficients
     """
    return sum(factorial(l)*binom(l,i)*\
               sum(Q[j]/sqrt(factorial(j[0])*factorial(j[1])) \
                   for j in const_ij[2*l]) for l in np.arange(i,m+1))


def write_SDPA_pb_upper(m,a,filename):
    """
    Aim
    ---
    Producing a text file for the SDP problem of computing an upper bound at 
    rank m in the hierarchy associated to the witness sum_k a_k |k><k|

    Parameters
    ----------
    m : integer
        rank in the hierarchy of SDP.
    a : list
        description of the witness.

    Returns
    -------
    Nothing
    But create a text file in the format ".dat-s" that can be later read by SDPA-GMP

    """          
    
    # Create list of indices (i,j) for which i+j= const 
    const_ij = [None]*(2*m+1)
    for i in range(m+1):
        for j in range(m+1):
            if const_ij[i+j] == None:
                const_ij[i+j] = []
            const_ij[i+j].append((i,j))

    # Create a picos problem
    D = pc.Problem()

    # Variables
    y = pc.RealVariable("y",1)
    mu = pc.RealVariable("mu",m+1)
    Q = pc.SymmetricVariable("Q",(m+1,m+1))
       
    # Objective
    D.set_objective("min",y)
    
    # Constraints
    for i in range(m+1):
        if i < len(a):
            D.add_constraint(y >= a[i]+mu[i])
            #D.add_constraint(mu[i] >= -1)
        else:
            D.add_constraint(y >= mu[i])
            #D.add_constraint(mu[i] >= -1)
        
    
    for i in range(m+1):
        D.add_constraint(get_Qi(Q,i,const_ij,m) == mu[i])
    
    D.add_constraint(Q >> 0)
    print(D)
    
    D.write_to_file(filename)
    

    
    
if __name__ == "__main__":
    # Location of the solver SDPA-GMP
    path_to_sdpa = '/home/pemeriau/Softwares/sdpa/sdpa-gmp-7.1.3/'
    
    if not os.path.isfile('param.sdpa'):
        os.system('cp '+path_to_sdpa+"param.sdpa"+" .")
    
    # Parameters
    t = 4 # rescaling coefficients. Won't affect the objective but might play an 
          # imporant role for the stability

    a = [1]  # vector defining the witness sum_k a_k |k><k| for which
                 # we will compute an upper bound an the threshold value.
                 # Careful, contrary to the conventionns of the main text, 
                 # vector a is written from Fock state 0 to Fock state len(a).
                 # Note that the level of the hierarchy should be higher than len(a).
                 # 
                 # Forall k, 0 <= a_k <= 1
                 # Normlisation: max_k a_k = 1
    
    m = 2 # rank in the hierarchy
    
    filename = "upper"
    for j in range(len(a)):
        filename += str(a[j])+"_"
    filename += "r"+str(m)+".dat-s"
    
    write_SDPA_pb_upper(m,a,filename)
    
    os.system(path_to_sdpa+"sdpa_gmp "+filename+" "+filename[0:-6]+".out")
    
