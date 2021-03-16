--- SA folder ---
This folder contains a couple of tools to perform sensitivity analysis.

In particular:
> SEIRD_PAWN.m: performs the PAWN (Pianosi-Wagener) method to evaluate the 
                dependence of the output upon the params, exploiting
                the KS norm
> SEIRD_PAWN_ic.m: performs the PAWN (Pianosi-Wagener) method to evaluate the 
                dependence of the output upon the params, exploiting
                the KS norm. Moreover, initial conditions can be passed
                as random variable
> SEIRD_VBSA.m: computes the Sobol' indices through Saltelli method.
                safe_R1.1 library is required.

Find safe_R1.1 at: https://www.safetoolbox.info
and save it in this current folder SA