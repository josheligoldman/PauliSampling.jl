using Printf

macro time_line(expr)
    quote
        local t = @elapsed local result = $(esc(expr))
        @printf("%-50s  %10.6f seconds\n", string($(QuoteNode(expr))), t)
        result
    end
end

macro alloc_line(expr)
    quote
        local a = @allocated local result = $(esc(expr))
        @printf("%-50s  %10d bytes\n", string($(QuoteNode(expr))), a)
        result
    end
end

function filter_no_XY(psum)
    nq = psum.nqubits
    filtered_terms = Base.filter(pair -> !containsXorY(pair.first), psum.terms)
    return PauliSum(nq, filtered_terms)
end

function get_one_inds!(one_inds, num_qinds, i)
    for j in 1:num_qinds
        if ((i >> (j-1)) & 1) == 1
            push!(one_inds, j)
        end
    end
    return
end

function apply_projector(psum, qind, x_k::Bool=false)
    # TODO: Make this a gate. 
    nq = psum.nqubits
    new_psum = PauliSum(ComplexF64, nq)
    sgn = (-1)^(x_k)
    B = PauliString(nq, :Z, qind, sgn)
    for (term, coeff) in psum
        pauli = getpauli(term, qind)
        if ispauli(pauli, 1) | ispauli(pauli, 2) # X or Y
            continue # discard since |x><x|X|x><x| = 0 and |x><x|Y|x><x| = 0 for x in {0, 1}
        end
        pstr = PauliString(nq, term, complex(coeff))
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

function qinds_in_order(i, num, nq)
    return i:min(i + (num-1), nq)
end

function qinds_random(i, num, nq)
    rng = i:nq
    return k ≤ length(rng) ? rand(rng, num; replace=false) : collect(rng)
end

qinds_base_func = (i, nq) -> qinds_in_order(i, 1, nq)

function random_psum(nq::Integer)
    labels = [:I, :X, :Y, :Z]
    total = 4^nq
    pstrs = Vector{PauliString}(undef, total)

    for idx in 0:total-1
        paulis = Vector{Symbol}(undef, nq)
        tmp = idx
        for i in 1:nq
            paulis[i] = labels[mod1(tmp, 4)]
            tmp ÷= 4
        end
        coeff = randn() + randn()*im
        pstrs[idx + 1] = PauliString(nq, paulis, 1:nq, coeff)
    end

    return PauliSum(pstrs)
end