--- SEIRD code ---
This library have been created to perform a complete workflow for the analysis 
of SEIRD model for epidemics.

In particular, the following steps can be performed:
- FwdUQ (Monte Carlo)
- Sensitivity analysis (PAWN, VBSA)
- InvUQ (MCMC-DRAM)

Each subfolder could also be run independently, paying attention to
define all the needed variables.
Moreover, a more specific _README.txt file is contained in each subfolder.
 
Our project was based on the paper: 
    C. Piazzola, L. Tamellini, R. Tempone, 
    "A note on tools for predictionunder uncertainty and identifiability
    of SIR-like dynamical systems forepidemiology"
 
suggesting us some aspects of the workflow and the data of the parameters.

WARNING: when calling MCMCstat, a problem with the function mad could occur.
        To overcome this issue, pay attention to include MCMCstat directory
        when performing inverse UQ and to exclude it when finished.
        All this is already implemented in example_main.m.
 
Needed packages:
- MCMCstat: https://mjlaine.github.io/mcmcstat/
- safe_R1.1: https://www.safetoolbox.info

M.Cesaratto, G.Pigani - march, 2021
 