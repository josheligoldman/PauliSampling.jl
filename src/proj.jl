function project(psum::PauliSum{TT, CT}, qind, bit::Bool) where {TT, CT}
    # TODO: Make this a gate. 
    nq = psum.nqubits
    new_psum = PauliSum(CT, nq)
    sgn = 1 - 2 * bit
    B = PauliString(nq, :Z, qind, sgn)
    for (term, coeff) in psum
        pauli = getpauli(term, qind)
        if ispauli(pauli, 1) | ispauli(pauli, 2) # X or Y
            continue # discard since |x><x|X|x><x| = 0 and |x><x|Y|x><x| = 0 for x in {0, 1}
        end
        pstr = PauliString(nq, term, coeff)
        # |x><x| p |x><x| = (1/2)(I + B) * p * (1/2)(I + B) = (1/4) * (p + pB + Bp + BpB)
        # Can do better, since anything that isn't an I among p, pB, Bp, BpB can be discarded as 
        # it will always be 0 in future trace calculations. 
        add!(new_psum, pstr)
        add!(new_psum, pstr * B)
        add!(new_psum, B * pstr)
        add!(new_psum, B * pstr * B)
    end
    mult!(new_psum, 1/4)
    return new_psum
end

function trace_project(psum::PauliSum{TT, CT}, qind, bit::Bool) where {TT, CT}
    # TODO: Make this non-allocating.
    nq = psum.nqubits
    aux_psum = PauliSum(CT, nq)
    s = 1 - 2 * bit
    for (pstr, coeff) in psum
        pauli = getpauli(pstr, qind)
        if ispauli(pauli, 1) || ispauli(pauli, 2) # X or Y
            continue
        end

        pstr = setpauli(pstr, 0, qind)
        new_coeff = coeff / 2
        if ispauli(pauli, 3)
            new_coeff *= s
        end

        add!(aux_psum, pstr, new_coeff)
    end
    return aux_psum
end
