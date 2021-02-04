#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementation of the hierarchy of UPPER bounds

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
m = n+20    #which rank in the hierarchy
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
D = pc.Problem()
D.options.solver = "mosek"  #change to any solver recognised by Picos (cvxopt for instance)

## Define variables
y = pc.RealVariable("y",1)
mu = pc.RealVariable("mu",m+1)
Q = pc.SymmetricVariable("Q",(m+1,m+1))
    

## Set objective
D.set_objective("min",y)

## Set constraints

# F_k >= 0
for i in range(m+1):
    if i==n:
        D.add_constraint(1+mu[i] <= y)
    else:
        D.add_constraint(mu[i] <= y)
    

# Sum_k (-1)^k mu_k L_k(tX) = X^TQX where X = (1,X/sqrt(1!),X^2/sqrt(2!),...)
for l in range(m+1):
    D.add_constraint(sum(Q[j]/sqrt((factorial(j[0])*factorial(j[1]))) for j in const_ij[2*l]) 
                     == sum( (-1)**(k+l)*binom(k,l)/factorial(l)*mu[k]*t**l   for k in np.arange(l,m+1)) )
    
# SDP constraint
D.add_constraint(Q >> 0)

## Save SDP as dat-s to run it later with SDPA
filename = 'upper_'+str(n)+'_'+str(m)+'dat-s'
D.write_to_file(filename)


## Solve SDP and print solution
solution = D.solve()

print(D.value)
print('----------------')
for i in range(len(mu.value)):
    print("F_",i,"    : ",round(mu.value[i],6))
