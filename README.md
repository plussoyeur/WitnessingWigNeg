# [Zeilberger's algorithm](Zeilberger_implemenation_feasibility.nb)

To run this notebook, you need to install RISC's algorithms for Mathematica (it will contain the program fastZeil of interest). Please find download and installation instruction [here](https://www3.risc.jku.at/research/combinat/software/ergosum/installation.html#download). You will need to contact [Peter Paule](https://risc.jku.at/m/peter-paule/) to be able to download the package. 

# Hierarchy of semidefinite programs ([upper](upperSDP.py) and [lower](lowerSDP.py))

These python codes were written using the interface provided by [Picos](https://picos-api.gitlab.io/picos/) which can be easily installed. You will further need a solver. A first option is to solve them directly by specifying a solver to Picos. You can install for instance [CVXOPT](https://cvxopt.org/) or [Mosek](https://www.mosek.com/) ([academic licenses](https://www.mosek.com/products/academic-licenses/) are available). A second, more precise option which allows to go higher in the hierarchy of SDP is to save the program with a 'dat-s' extension and then solve the program using [SDPA-GMP](https://sourceforge.net/projects/sdpa/files/sdpa-gmp/) which is an SDP solver with precision arithmetic.

Feel free to contact any of the authors should you have any difficulties in implementing the codes.
