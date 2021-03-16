--- InvUQ folder ---
This folder aims to perform inverse analysis on the SEIRD model

In particular:
> SERIDlsq.m: implements a Least Square estimator, calling SEIRDss.m to
                perform error minimization
> SEIRDsynt.m: generates a set of syntethic data for the specified compartements
                and with the desired variance, centered on the results
                obtained by the "true" parameters
> SEIRDmcmc_ic: performs teh MCMC (DRAM) algorithm itself, relying on
                the library MCMCstat library

Find MCMCstat at: https://mjlaine.github.io/mcmcstat/
and save it in this current folder (InvUQ)