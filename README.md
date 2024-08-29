# xStressorsStabBEFW

This repo contains the simulation, the code, analysis and the manuscript of the
study: "Danet, A, Kéfi, S, Johnson, T, F & Beckerman, A, P. Response diversity is a major driver of temporal stability in
complex food webs."

This work relies on simulations the dynamics of complex food-webs using the
bioenergetic food-web model.

The bioenergetic model is implemented in the
[EcologicalNetworksDynamics.jl](https://github.com/BecksLab/EcologicalNetworksDynamics.jl)
julia package and presented in [Lajaaiti et al. (2024)](https://www.biorxiv.org/content/10.1101/2024.03.20.585899v1).

For this study, we developped a stochastic extension of this model by allowing
species mortality rates to have a stochastic component. This extension has been
developped in a fork of the original package [in my
repo](https://github.com/alaindanet/BEFWM2/tree/vasseur_fox).

This model is basically an extension of the model of [Vasseur & Fox
(2007)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1461-0248.2007.01099.x)
to complex food-webs and with all species having stochastic mortality rates.

The simulations were thus run in `Julia`, but the statistical analysis were done
in `R`.

# The steps and organisation

## The simulations

All simulations scripts lives in the `scripts` folder. The parameters of the
simulations and the generation of food-web topology is described in
`param_comb_zc.jl`

We ran simulations with this script: `simulation_args.jl`. As it indicates, this script
read the arguments provided at execution of the script. This methods allowed to
easily run sensibility analysis while keeping the same master script.

We ran 4 types of simulations, presented in Fig. S1 of the manuscript:
- Default simulations
- Simulations without standardising the carrying capacity of the producers, i.e.  
without dividing carrying capacity by the number of primary producers
- Simulations with rebuilding the preferences of the consumers after removing
  disconnected/dead species
- Simulations with fixed death rates instead of allometric

The timeseries of the simulations were used to compute stability and food-web
metrics in a separated repository using R: [alaindanet/pre_process_xStressorsStabBEFW](https://github.com/alaindanet/pre_process_xStressorsStabBEFW)

## Analysis and manuscript

The statistical analysis, the figures, and the tables were done in this repo
using the [targets](https://books.ropensci.org/targets/) R package, which is
a pipeline tool ensuring that the results are up-to-date and reproducible. The
pipeline is defined in `_target.Rmd`.

The manuscript and supplementary material are written in Rmardown, using
`targets` to print the numbers and the figures. The are located in the `paper`
folder.


