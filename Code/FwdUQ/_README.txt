--- FwdUQ folder ---
It can be used to perform a (crude) Monte Carlo analysis of the SEIRD model.

In particular:
> SEIRDmc.m: it's intended to receive fixed initial conditions
> SEIRDmc_ic.m: it's intended to receive the bounds between which initial
                conditions are included (uniformely distributed)

Both of them perform nMC evaluations of the SEIRD model, dispaying (optionally)
the dynamics and the PDFs distributions of QoIs