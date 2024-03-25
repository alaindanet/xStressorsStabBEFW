println("Rundir is $(pwd())")

import Pkg
using Distributed, Serialization, ArgParse

@everywhere include("../src/misc.jl")

args = parse_commandline()

first_sim = args["first_sim"]
last_sim = args["last_sim"]

if args["d"] == "nothing" 
    args["d"] = nothing
else
    args["d"] = parse(Float64, args["d"])
end

print_argument()

#ncpu = maximum([length(Sys.cpu_info()), 15])
ncpu = 20

#Dir
proj_dir = expanduser("~/xStressorsStabBEFW")

#Flag enables all the workers to start on the project of the current dir
flag = "--project=" * proj_dir
println("Workers run with flag: $(flag)")
addprocs(ncpu - 1, exeflags=flag)
#addprocs(5, exeflags=flag)
println("Using $(ncpu -1) cores")


@everywhere import Pkg, Random.seed!
@everywhere using DifferentialEquations, EcologicalNetworksDynamics, SparseArrays
@everywhere using LinearAlgebra, DataFrames
@everywhere using Distributions, ProgressMeter
@everywhere using StatsBase, CSV, Arrow
@everywhere include("../src/stochastic_mortality_model.jl")
@everywhere include("../src/sim.jl")
@everywhere include("../src/interaction_strength.jl")
