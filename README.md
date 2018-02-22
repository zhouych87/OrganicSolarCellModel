# Organic Solar Cell Model
Towards predicting the power conversion efficiencies of organic solar cells from donor and acceptor molecule structures
For the detials please check my paper:http://pubs.rsc.org/en/content/articlelanding/2014/TC/C7TC05290A#!divAbstract
Please cite my paper if you use my codes. 

# Codes
* gfilm (split PBC thin film into molecules and dimers; output molecule connection information; output Gaussian input files to calculate electronic coupling)
* calv (read Gaussian ouput log file and calculate orbital couplings)
* rw (random walk simulation)
* solar.out (solar cell numerical simulation)
* g09log (read HOMO and LUMO number,other usage)

# How to use it?
In my paper, Figure 2 shows the overall processes of this multi-scale simulation. Given donor and acceptor molecular structures, we calculate the reorganization energies and free energies by GAUSSIAN 09 package. We can also abstract force field parameters from this calculation. These force field parameters will be used in the mesoscopic molecule dynamics (MD) simulations. Second, we insert donor and acceptor according to experiment ratios into a 10 nm $\times$ 10 nm $\times$ 10 nm box. Using calculated force field to perform MD simulations at 400 K and then cool it the 300 K for 10 ns and then continues another 10 ns after reached equilibrium. These MD simulations are carried out by GROMACS. Third, using the home-made code 'gfilm' to split thin film models into molecules and dimers, which are written with the Gaussian input format. At the same time, geometric center of the molecules, connection information and a script named with 'run' to call Gaussian and 'calv' are generated. Gaussian calculations output HOMO, LUMO and molecular orbitals. Then'calv' reads the Gaussian output log file and calculates electronic couplings. These calculated couplings together with the calculated free energies and reorganization energies are implemented into Marcus theory to estimate the charge transfer (CT) rate between every dimer. Using calculated CT rates and thin film geometry information, we perform random walk (RW) simulations, in which we obtain electron and hole diffusion coefficients. Corresponding mobilities are calculated by the Einstein relation $\mu=\frac{eD}{kT}$. Inputting these calculated diffusion coefficient and mobilities of electron and hole, the band gap, effective density states of the 'conduction band' and 'valence band', and three experiment parameters ($\alpha$, $\tau$ and $\epsilon$) into the numerical model, we obtain the J-V curve of the solar cell. PCEs, FF, V$_{oc}$ and J$_{sc}$ are read from the J-V curve. 
