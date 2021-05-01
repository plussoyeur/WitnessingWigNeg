# [Zeilberger's algorithm](https://github.com/plussoyeur/WitnessingWigNeg/blob/main/Zeilberger_implemenation_feasibility.nb)

To run this notebook, you need to install RISC's algorithms for Mathematica (it will contain the program fastZeil of interest). Please find download and installation instruction [here](https://www3.risc.jku.at/research/combinat/software/ergosum/installation.html#download). You will need to contact [Peter Paule](https://risc.jku.at/m/peter-paule/) to be able to download the package. 

# Hierarchy of semidefinite programs ([upper](https://github.com/plussoyeur/WitnessingWigNeg/tree/main/upper) and [lower](https://github.com/plussoyeur/WitnessingWigNeg/tree/main/lower))
These python codes were written using the interface provided by [Picos](https://picos-api.gitlab.io/picos/) which can be easily installed. You will further need a solver. A first option is to solve them directly by specifying a solver to Picos. You can install for instance [CVXOPT](https://cvxopt.org/) or [Mosek](https://www.mosek.com/) ([academic licenses](https://www.mosek.com/products/academic-licenses/) are available). A second, more precise option which allows to go higher in the hierarchy of SDP is to save the program with a 'dat-s' extension and then solve it using [SDPA-GMP](https://sourceforge.net/projects/sdpa/files/sdpa-gmp/) which is an SDP solver with precision arithmetic. Simple examples with any solver can be ran for the upper hierarchy [here](https://github.com/plussoyeur/WitnessingWigNeg/blob/main/upper/upperSDP.py) and for the lower hierarchy [here](https://github.com/plussoyeur/WitnessingWigNeg/blob/main/lower/lowerSDP.py). 
Codes here are implemented with SDPA-GMP. Files of the form `gen_SDPA` ([upper](https://github.com/plussoyeur/WitnessingWigNeg/blob/main/upper/gen_SDPA_upper.py) and [lower](https://github.com/plussoyeur/WitnessingWigNeg/blob/main/lower/gen_SDPA_lower.py)) exploits Picos to generate a file with extension `dat-s` readable by SDPA-GMP. Files `solve_` ([upper](https://github.com/plussoyeur/WitnessingWigNeg/blob/main/upper/solve_upper.py) and [lower](https://github.com/plussoyeur/WitnessingWigNeg/blob/main/lower/solve_lower.py)) runs the entire hierarchy until SDPA-GMP raises a flag of non-optimality. 

A list of 331 precomputed lower bounds and upper bounds on the threshold values for witnesses with n=3 is given [here](https://github.com/plussoyeur/WitnessingWigNeg/blob/main/thresholdvalues_a1a2a3.txt).

Feel free to contact any of the authors should you have any difficulties in implementing the codes.
