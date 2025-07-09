function sample_bitstring(psum; prob_method::Symbol = :approx)
    nq = psum.nqubits
    running_psum = psum
    bitstring = BitVector(undef, nq)

    if prob_method == :coset
        f = coset_normalizing_factor(running_psum)
    end

    for qind in 1:nq
        if prob_method == :coset
            proj0 = apply_projector(running_psum, qind, false)
            bitstring[qind] = false
            f0 = coset_normalizing_factor(proj0; num_set_bits=qind, set_bits=@view bitstring[1:qind])
            p0 = f0 / f

            if rand() <= p0
                bitstring[qind] = false
                f = f0
                running_psum = proj0
            else
                bitstring[qind] = true
                f = f - f0
                running_psum = running_psum - proj0
            end

        elseif prob_method == :approx

            p0 = normalized_prob(running_psum, qind, 0)        
            
            if rand() <= p0
                bitstring[qind] = false
            else
                bitstring[qind] = true
            end
            running_psum = apply_projector(running_psum, qind, bitstring[qind])

        else
            error("Unknown prob_method: $prob_method. Use :coset or :approx.")
        end
    end

    return bitstring
end