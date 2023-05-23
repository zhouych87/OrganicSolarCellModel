
Intermolecular 3D-MoRSE descriptors for fast and accurate prediction of electronic couplings in organic semiconductors
T

# Codes
* gfilm (split PBC thin film into molecules and dimers; output molecule connection information; output Gaussian input files to calculate electronic coupling)
* calv (read Gaussian ouput log file and calculate orbital couplings)
* rw (random walk simulation to calculation hole and electron mobility of OSC)
* solar.out (solar cell numerical simulation)
* g09log (give HOMO and LUMO levels, other usage)
* g16log.f90 compile by "gfortran g16log.f90 -o g16log". Use as "g19log xxxx.log 1", will output molecular orbital levels in the range of -7.0 eV -2.0 eV 

# How to use it?
In my paper, Figure 2 shows the overall processes of this multi-scale simulation. Given donor and acceptor molecular structures, we calculate the reorganization energies and free energies by GAUSSIAN 09 package. We can also abstract force field parameters from this calculation. These force field parameters will be used in the mesoscopic molecule dynamics (MD) simulations. Second, we insert donor and acceptor according to experiment ratios into a 10 nm $\times$ 10 nm $\times$ 10 nm box. Using calculated force field to perform MD simulations at 400 K and then cool it the 300 K for 10 ns and then continues another 10 ns after reached equilibrium. These MD simulations are carried out by GROMACS. Third, using the home-made code 'gfilm' to split thin film models into molecules and dimers, which are written with the Gaussian input format. At the same time, geometric center of the molecules, connection information and a script named with 'run' to call Gaussian and 'calv' are generated. Gaussian calculations output HOMO, LUMO and molecular orbitals. Then'calv' reads the Gaussian output log file and calculates electronic couplings. These calculated couplings together with the calculated free energies and reorganization energies are implemented into Marcus theory to estimate the charge transfer (CT) rate between every dimer. Using calculated CT rates and thin film geometry information, we perform random walk (RW) simulations, in which we obtain electron and hole diffusion coefficients. Corresponding mobilities are calculated by the Einstein relation $\mu=\frac{eD}{kT}$. Inputting these calculated diffusion coefficient and mobilities of electron and hole, the band gap, effective density states of the 'conduction band' and 'valence band', and three experiment parameters ($\alpha$, $\tau$ and $\epsilon$) into the numerical model, we obtain the J-V curve of the solar cell. PCEs, FF, V$_{oc}$ and J$_{sc}$ are read from the J-V curve. 

* Build donor and accepter molecules strucutres and calculate force field and reorganization energies. (Manually)
* Perform MD simulation with GROMACS
* Use gfilm to split the final '.gro' file, it will generate: connection.dat,grp,tcon,run,*.com 
* run the 'run' file to call Gaussian, please keep in mind that this may take very very very long time due to large amount of molecules and dimers!!!! Obitals are calculated by Gaussian, calv will read these orbital and calculate orbital coupling of dimers. Gaussian output file "*.log" will be deleted because each log file are usually several Gb, which will use up your disk.     
* run the 'grp' file to call grep to collect coupling data
* run rw to obtain the diffusion coefficient and mobility, please edit the j.out file to implement reorganization energies
* input parameter into 'in', run solar.out  



# Please edit the 'in' file in each folder
'in' is the input parameter for the program. Please revise it according to your data.

# Please let me know if there was any bug. Thank you! 
My email: zhouych87@gmail.com.
