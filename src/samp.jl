function sample_bitstring(nq, psum, marginal_func = bayes_marginal)
    running_psum = psum
    bitstring = BitVector(undef, nq)
    for qind in 1:nq
        p0 = marginal_func(running_psum, qind, false)
        bitstring[qind] = rand() <= p0
        running_psum = project(running_psum, qind, bitstring[qind])
    end
    return bitstring
end