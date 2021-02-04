#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: pemeriau

Implementation of the hierarchy of LOWER bounds

Our targeted witness is the fidelity with a single Fock state |n> but the code
can be easily extended to the case of linear combination with several Fock states
"""

## import
import picos as pc
import numpy as np
from math import factorial, sqrt
from scipy.special import binom


## Parameters
n = 3       #specify the witness |n><n|
m = n+2   #which rank in the hierarchy
t = 4       #dilation in the Wigner function: \sum_k (-1)^k F_k L_k(tX)


## Useful for later
# Create list of indices (i,j) for which i+j = constant
const_ij = [None]*(2*m+1)
for i in range(m+1):
    for j in np.arange(m+1):
        if const_ij[i+j] == None:
            const_ij[i+j] = []
        const_ij[i+j].append((i,j))



## Define problem
P = pc.Problem()
P.options.solver = "mosek"  #change to any solver recognised by Picos (cvxopt for instance)

## Define variables
F = pc.RealVariable("F",m+1)
Q = pc.SymmetricVariable("Q",(m+1,m+1))

## Set objective
P.set_objective("max",F[n])

## Set constraints

# F_k >= 0
for i in range(m+1):
    P.add_constraint(F[i] >= 0)
    
# Sum_k F_k = 1
P.add_constraint(sum(F[i] for i in range(m+1)) == 1)

# Sum_k (-1)^k F_k L_k(tX) = X^TQX where X = (1,X/sqrt(1!),X^2/sqrt(2!),...)
for l in range(m+1):
    P.add_constraint(sum(Q[j]/sqrt((factorial(j[0])*factorial(j[1]))) for j in const_ij[2*l]) 
                     == sum( (-1)**(k+l)*binom(k,l)/factorial(l)*F[k]*t**l   for k in np.arange(l,m+1)) )
for l in np.arange(1,m+1):
    P.add_constraint(sum(Q[j]/sqrt((factorial(j[0])*factorial(j[1]))) for j in const_ij[2*l-1]) 
                     == 0 )
    
# SDP constraint
P.add_constraint(Q >> 0)

## Save SDP as dat-s to run it later with SDPA
filename = 'lower_'+str(n)+'_'+str(m)+'dat-s'
P.write_to_file(filename)


## Solve SDP and print solution
solution = P.solve()

print(P.value)
print('----------------')
for i in range(len(F.value)):
    print("F_",i,"    : ",round(F.value[i],6))
