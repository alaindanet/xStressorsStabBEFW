#
mkdir("src")

cp("/home/alain/Documents/post-these/sheffield/befwm2_dev/Manifest.toml", "/home/alain/Documents/post-these/sheffield/xStressorsStabBEFW/Manifest.toml")
cp("/home/alain/Documents/post-these/sheffield/befwm2_dev/Project.toml", "/home/alain/Documents/post-these/sheffield/xStressorsStabBEFW/Project.toml")

mkdir("use_case")
touch("src/minmax.jl")
touch("draft.jl")

Pkg.add("Plots")

using Pkg
Pkg.dev("~/Documents/post-these/sheffield/BEFWM2_fork/")
