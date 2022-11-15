#
mkdir("src")

cp("/home/alain/Documents/post-these/sheffield/befwm2_dev/Manifest.toml", "/home/alain/Documents/post-these/sheffield/xStressorsStabBEFW/Manifest.toml")
cp("/home/alain/Documents/post-these/sheffield/befwm2_dev/Project.toml", "/home/alain/Documents/post-these/sheffield/xStressorsStabBEFW/Project.toml")

mkdir("use_case")
touch("src/minmax.jl")
touch("draft.jl")

Pkg.add("Plots")

using Pkg
Pkg.dev("../BEFWM2_fork/")

mkdir("res")

mkdir("report")
touch("report/project_definition.Rmd")

symlink("/home/alain/Documents/post-these/references.bib", "report/references.bib")

