<h1>GPshared</h1>

<p>This project is a c/c++ shared library and python wrapper to be used for Gaussian Process Global Optimisation. The GPpython folder contains a separate python project for Bayesian Optimization</p>
<p>The cpp library and python wrapper implement Gaussian processes with observations and inference up to second derivatives in multiple axes, and posterior inference over sets of equally weighted hyperparameter draws</p>
<h3>Usage:</h3>
<p>running main.py should run some demos and show the plots in turn. The demos are
<UL>
    <LI> making some draws from a GP</LI>
    <LI> showing observation and inference of derivatives</LI>
    <LI> drawing hyperparameter samples from from the posterior</LI>
    <LI> optimizing branin with zero mean noise but with a cost relation on variance, variable fidelity outperforms fived for PES and EI</LI>
    <LI> optimizing a synthetic fn with biased results but low noise and a cost relation on how close to true the result is. PES with env variable outperforms fixed at any point.</LI>
</UL></p>

<h3>Dependencies (that I can think of):</h3>
<p>numpy/scipy matplotlib pathos.multiprocessing(shouldn't be needed for demos as the result files are in the cache) DIRECT tqdm</p>
    
<h3>Files:</h3>
<p>
    In c lib
<UL>
    <LI> libGP.cpp: main file for cpp library, maintains a vector of GPs and handles calls on them, also a function to return log likelihood only</LI>
    <LI> GPsimple.cpp: GP implementation, methods for infering mean only, diagonal variance or full covariance</LI>
    <LI> direct.cpp: direct from external source</LI>
    <LI> bayesutils.cpp: EI, LCB and logEI functions</LI>
    <LI> hypsearch.cpp: fn to run direct to find MLE or MAP hyperparameters given data (and prior for MAP)</LI>
    <LI> kernel.cpp: has or imports all the kernel functions, taking numnbers correspoinding to derivatives in each axis of each x1 and x2, converters for hyperparameters from  natural (length) to hte form used in kernel (1/l**2) and to log space for searching</LI>
    <LI> matern.cpp: the matern 5.2 kernel</LI>
    <LI> misctool.cpp: normal pdf and cdf, EI and logEI draw samples from a covariance matrix or cholesky decomp.</LI>
    <LI> simplekernels.cpp: some sums and productts of other kernels, derivatives not implemented.</LI>
   </UL>
    In python lib
<UL>
    <LI>   GPdc.py: imports the c library and provides interface</LI>
    <LI>   slice.py: slice samples adapted from 3rd party</LI>
    <LI>   search.py: wrapper for search algorithms, sets up parameters and runs the search</LI>
    <LI>   eprop.py: expectation propagation for hard and soft inequalities on multivatiate normals</LI>
    <LI>   optutils.py: selection of synthetic objectives and implementation of all the searches, mostly inheriting from a base</LI>
    <LI>   ESutils.py: draws support points and mins from support given a GP, draws hyperparameters from data given prior</LI>
    <LI>   PES.py: implementation of PESaquisition function, regular and env variable versions.</LI>
    <LI>   test* : tests for all the above</LI>
    <LI>   exp*  : experiments according to comments in header</LI>
    </UL>
    
    </p>

<h3>Using External Code:</h3>
<p>DIRECT implementation<BR> Jones, D.R., Perttunen, C.D., Stuckman, B.E.: Lipschitzian optimization without 
!    the Lipschitz constant. J. Optim. Theory Appl. 79(1), 157–181 (1993)</p>