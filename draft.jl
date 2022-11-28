using Revise
using BEFWM2
using DifferentialEquations
using Statistics
using DataFrames
using Plots
using Debugger
include("src/minmax.jl")

Iₖ = ωₖ * x * y / B₀


# Prepare your names (Symbols rather than Strings)
names = (:a, :b, :c)
# Prepare your values (no need to `collect` them yet).
as = 1:2
bs = 3:4
cs = 5:6
# Iterate on cartesian product (the only necessary `collect`ion happens here).
for product in Iterators.product(as, bs, cs)

    # Here is one element of your product.
    println("\nproduct: $product {$(typeof(product))}")

    # Construct a dict with (Symbol => value) pairs.
    dict = Dict(key => value for (key, value) in zip(names, product))
    println("dict: $dict {$(typeof(dict))}")

    # "Splat" the dict into a named tuple.
    named_tuple = (;dict...)
    println("named_tuple: $named_tuple {$(typeof(named_tuple))}")
end

# One-liner:
map(p -> (;Dict(k => v for (k, v) in zip(names, p))...), Iterators.product(as, bs, cs))[:]

using BEFWM2
BEFWM2.PositiveDomain()
