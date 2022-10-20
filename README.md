# Code for the article: On the equation x²+dy⁶=z^p for square-free 1≤d≤20.

This git repository contains files related to the computational part of the article: On the equation x²+dy⁶=z^p for square-free 1≤d≤20. It contains the codes used to run "Mazur's trick", ir order to discard the newforms, and also the code to compute Ellenberg's bounds. It also includes the codes to compute Table 2.1 and the scripts to get with the minimal models of the curves. You can find the outputs in "Outputs.txt".

# Contents:

**• mazur.m** : To run Mazur's trick, in Magma.

**• Multi-Frey.m** : Multi-Frey code to the case d=4,  in Magma.

**• ApCoefficients.m** : To compute coefficients of newforms, in Magma.

**• Ellenberg.gp** : To compute Ellenberg's bounds for Proposition 3.2, in PARI/GP.

**• Cased4.m** :  The code used in Theorem A.8, in Magma.

**• Cased5.m** :  The code used in Theorem 3.3, in Magma.

**• Cased7.m** :  The code used in Theorem 3.4, in Magma.

**• Cased11.m** :  The code used in Theorem 3.5, in Magma.

**• Cased13.m** :  The code used in Theorem 3.6, in Magma.

**• Cased15.m** :  The code used in Theorem 3.7, in Magma.

**• Cased19.m** :  The code used in Theorem 3.8, in Magma.

**• Cased19.gp** :  The code used in Theorem 3.8, in PARI/GP.

**• coef5.txt** :  A file to run Cased19.gp (an output of Cased19.m).

**• coef11.txt** :  A file to run Cased19.gp (an output of Cased19.m).

**• pol5.txt** : A file to run Cased19.gp (an output of Cased19.m).

**• Table.m** :  The code used to compute Table 2.1, in Magma.

**• Minimal_model** :  This folder contains the files "Curves_d=n.sage" for n=5,7,11,15., with the minimal models of the curves that appear in the paper (in Sage).

**• Minimal_model** :  This folder contains the files "Multi-Frey.m" and "Codesd4.m", needed to run the file "Cased4.m", in Magma.

**• Cased4.m** : The code used in Theorem A.8, in Magma.

**• Outputs.txt** :  A file that shows the outputs of "Casedn.m/gp", for n=5,7,11,13,15,19. 
