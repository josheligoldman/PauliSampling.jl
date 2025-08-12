function sample_bitstring(psum, marginal_func = bayes_marginal, proj_func = trace_project)
    nq = psum.nqubits
    running_psum = psum
    bitstring = BitVector(undef, nq)
    for qind in 1:nq
        p0 = marginal_func(running_psum, qind, false)
        bitstring[qind] = rand() <= p0
        running_psum = proj_func(running_psum, qind, bitstring[qind])
    end
    return bitstring
end