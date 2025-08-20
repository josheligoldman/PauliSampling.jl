function sample_bitstring(psum, ::Type{UT}; marginal_func=bayes_marginal, proj_func=trace_project) where {UT <: Unsigned}
    nq = psum.nqubits
    running_psum = psum
    bitstring = zero(UT)
    for i in 1:nq
        p0 = marginal_func(running_psum, i, false)
        b_i = rand() ≥ p0
        bitstring = set_bit(bitstring, i, b_i)
        running_psum = proj_func(running_psum, i, b_i)
    end
    return bitstring
end