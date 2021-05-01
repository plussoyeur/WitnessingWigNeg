#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:25:15 2020

@author: pemeriau

In this file, we run the lower hierarchy until it stops having an optimility flag 
thus obtaining the highest numerical optimal value.
"""

# Import useful libraries. 
import re
import os
from gen_SDPA_lower import write_SDPA_pb_lower


# Location of the solver SDPA-GMP
path_to_sdpa = '/home/pemeriau/Softwares/sdpa/sdpa-gmp-7.1.3/'
    
if not os.path.isfile('param.sdpa'):
    os.system('cp '+path_to_sdpa+"param.sdpa"+" .")
  
    
def compute_hierarchy_lower(a):
    """
    Parameters
    ----------
    a : list
        Vector defining the witness sum_k a_k |k><k| for which
        It will compute a lower bound on the threshold value.
        Careful, contrary to the conventionns of the main text, 
            vector a is written from Fock state 0 to Fock state len(a).
        Note that the level of the hierarchy should be higher than len(a).
                  
        Forall k, 0 <= a_k <= 1
        Normlisation: max_k a_k = 1

    Returns
    -------
    Value of the higher rank in the hierarchy where optimality is certified
        together with input vector a.
    Last output file before computation is not optimal will be saved.
    """

    m = len(a) # start hierarchy at len(a)
    flag = "pdOPT"

    while flag == "pdOPT":
        # Create filename
        if m == len(a):
            filename = "lower"
            for j in range(len(a)):
                filename += str(a[j])+"_"
            filename += "r"+str(m)+".dat-s"
        
        if m != len(a):
            filename = "lower"
            filename_old = "lower"
            for j in range(len(a)):
                filename += str(a[j])+"_"
                filename_old += str(a[j])+"_"
            filename += "r"+str(m)+".dat-s"
            filename_old += "r"+str(m-1)+".dat-s"
            filename_old_out = filename_old[0:-6]+".out"
        
        filename_out = filename[0:-6]+".out"
    

         
        # Create SDP file readable by SDPA-GMP
        write_SDPA_pb_lower(m,a,filename)
        # Solve thanks to SDPA-GMP
        os.system(path_to_sdpa+"sdpa_gmp "+filename+" "+filename_out)
    
        with open(filename_out) as file:
            line = file.readline()
            while line:
                line=file.readline()
                if re.search("phase.value",line):
                    flag = line.split()[-1]
                if re.search("objValPrimal",line):
                    objP = float(line.split()[-1])

        # If optimality flag OK, delete old file otherwise delete current
        if m == len(a): #in case first rank already not optimality flag
            objP_stored = objP
        if m != len(a):
            if flag == "pdOPT":
                os.system('rm '+filename_old+' '+filename_old_out)
                objP_stored = objP
            else:
                os.system('rm '+filename+' '+filename_out)
            
        # Increment hierarchy
        m +=1
    
    return (a,objP_stored)

    
if __name__ == "__main__":
    a = [0,0,0,1]  # vector defining the witness sum_k a_k |k><k| for which
                 # we will compute an upper bound an the threshold value.
                 # Careful, contrary to the conventionns of the main text, 
                 # vector a is written from Fock state 0 to Fock state len(a).
                 # Note that the level of the hierarchy should be higher than len(a).
                 # 
                 # Forall k, 0 <= a_k <= 1
                 # Normlisation: max_k a_k = 1
    val = compute_hierarchy_lower(a)
