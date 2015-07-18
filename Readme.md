# SBFL
SBFL is a MATLAB implementation of the split Bregman method for fused Lasso. Please refer [a paper](http://www.sciencedirect.com/science/article/pii/S0167947310004093) (Split Bregman method for large scale fused Lasso, Ye et al., CSDA 2011) for details.<br />
Copyright (c) 2015 Taehoon Lee

# Functions
<li> SBFL.m: a linear system embedded in the SB-iteration is solved by matrix inversion. </li>
<li> SBFL_PCG.m: the liner system is solved by preconditioned conjugate gradients (PCG). </li>
<li> SBFL_CGLS.m: the liner system is solved by conjugate gradient for least squares problems (CGLS). </li>
<li> SBFL_PCGLS.m: the liner system is solved by preconditioned conjugate gradient for least squares problems (PCGLS). </li>

# Getting Started
An example is available on testSBFL.m.