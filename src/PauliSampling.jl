module PauliSampling

using PauliPropagation
using F2Algebra
using Hadamard
using Combinatorics
using BenchmarkTools
using Random
using LinearAlgebra
using LinearMaps
using Arpack
using Plots

include("utils.jl")
export 
    get_dist

include("probs.jl")
export 
    approximate_prob,
    projection_prob

include("dists.jl")

include("sampling.jl")
export 
    sample_bitstring

include("state.jl")
export 
    zero_state,
    build_circuit

# export coset_normalizing_factor


end # module PauliSampling
