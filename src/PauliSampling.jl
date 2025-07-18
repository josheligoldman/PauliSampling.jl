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

# Extend Base methods for PauliFreqTracker
import Base: real, imag, abs, complex, convert, float

real(p::PauliFreqTracker) = real(p.coeff)
imag(p::PauliFreqTracker) = imag(p.coeff)
abs(p::PauliFreqTracker) = abs(p.coeff)
float(p::PauliFreqTracker) = float(p.coeff)
complex(p::PauliFreqTracker) = complex(p.coeff)
convert(::Type{ComplexF64}, p::PauliFreqTracker) = convert(ComplexF64, p.coeff)

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

include("analysis.jl")
export
    tvd,
    kl_div,
    avg_metrics_for_weight,
    bit_marginals,
    shannon_entropy

end # module PauliSampling
