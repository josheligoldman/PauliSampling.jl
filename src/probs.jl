function normalizing_factor(psum::PauliSum{TT, CT}) where {TT, CT}
    nq = psum.nqubits
    acc = 0
    for one_inds in powerset(1:nq)
        acc += abs(overlapwithcomputational(psum, one_inds))
    end
    return acc
end

function faster_normalizing_factor(psum::PauliSum{TT, CT}) where {TT, CT}
    nq = psum.nqubits
    f_psum = filter_no_XY(psum)
    return normalizing_factor(f_psum)
end

function coset_normalizing_factor(psum::PauliSum{TT, CT}; to_filter=false, num_set_bits=0, set_bits=BitVector()) where {TT, CT}
    nq = psum.nqubits
    if to_filter
        f_psum = filter_no_XY(psum)
    end
    num_pstrs = length(f_psum)

    Z_matrix = BitMatrix(undef, num_pstrs, nq)
    coeffs = Vector{ComplexF64}(undef, num_pstrs)
    for (i, (term, coeff)) in enumerate(f_psum)
        for j in 1:nq
            pauli = getpauli(term, j)
            Z_matrix[i, j] = ispauli(pauli, 3)
        end
        coeffs[i] = coeff
    end

    @views Z_matrix_a = Z_matrix[:, 1:num_set_bits]
    Z_matrix_b = Z_matrix[:, num_set_bits+1:nq]

    y = f2_mul(Z_matrix_a, set_bits)

    _, rank, pivots = F2Algebra.rref(Z_matrix_b)
    @views basis = Z_matrix_b[:, pivots]

    acc = 0.0
    v = zero(Z_matrix[:, 1])
    for i in 0:(2^rank - 1)
        fill!(v, false)
        for j in 0:(rank - 1)
            if ((i >> j) & 1) == 1
                v .⊻= @view basis[:, j + 1]
            end
        end

        tmp = 0.0 + 0.0im
        for k in 1:num_pstrs
            sgn = ifelse(v[k] ⊻ y[k], -1.0, 1.0)
            tmp += coeffs[k] * sgn
        end
        acc += abs(tmp)
    end
    acc *= 2^(nq - rank - num_set_bits)
    return acc
end

function compute_grouped_traces(psum::PauliSum{TT, CT}, qinds::Vector{<:Integer}) where {TT, CT}
    # TODO: I don't think I'm normalizing correctly, but it's still giving correct results. Why?
    num_qinds = length(qinds)
    coeffs = Vector{ComplexF64}(undef, 2^num_qinds)
    one_inds = Int[]
    sizehint!(one_inds, num_qinds)
    Z_str = fill(:Z, num_qinds)

    for i in 0:(2^(num_qinds) - 1)
        empty!(one_inds)
        get_one_inds!(one_inds, num_qinds, i)
        @views Z_view = Z_str[1:length(one_inds)]
        @views qinds_view = qinds[one_inds]
        coeffs[i + 1] = getcoeff(psum, Z_view, qinds_view)
    end
    return ifwht(coeffs)
end

function grouped_normalizing_factor(psum::PauliSum{TT, CT}, qinds::Vector{Integer}) where {TT, CT}
    vals = compute_grouped_traces(psum, qinds)
    return sum(abs, vals)
end

function projection_prob(psum::PauliSum{TT, CT}, one_bit_inds::Vector{<:Integer}) where {TT, CT}
    numerator = abs(overlapwithcomputational(psum, one_bit_inds))
    denominator = coset_normalizing_factor(psum)

    return numerator / denominator
end

function projection_prob(psum::PauliSum{TT, CT}, x::BitVector) where {TT, CT}
    return projection_prob(psum, findall(x))
end

function projection_prob(psum::PauliSum{TT, CT}, x::Integer) where {TT, CT}
    nq = psum.nqubits
    return projection_prob(psum, BitVector(digits(x, base=2, pad=nq)))
end

function projection_marginal(psum::PauliSum{TT, CT}, idx::Integer, bit::Bool; to_filter=true, num_set_bits=0, set_bits=BitVector()) where {TT, CT}
    b = coset_normalizing_factor(psum; to_filter=to_filter, num_set_bits=num_set_bits, set_bits=set_bits)
    proj_psum = apply_projector(psum, idx, bit)
    a = coset_normalizing_factor(proj_psum; to_filter=to_filter, num_set_bits=num_set_bits, set_bits=set_bits) # TODO: Could add the set bit idx.

    return a / b
end

function grouped_marginal(psum::PauliSum{TT, CT}, idx::Integer, bit::Bool, qinds::Vector{<:Integer}) where {TT, CT}
    # idx must be the first element in qinds. Could get rid of this assumption, but probably not worth it. 
    @assert qinds[1] == idx
    vals = compute_grouped_traces(psum, qinds)
    p0_unnormalized = sum(abs, @view vals[1:end÷2])
    p1_unnormalized = sum(abs, @view vals[end÷2 + 1:end])
    return (bit == 0 ? p0_unnormalized : p1_unnormalized) / (p0_unnormalized + p1_unnormalized)
end

function grouped_prob(psum, x::BitVector; qinds_func=qinds_base_func)
    running_psum = psum
    nq = running_psum.nqubits
    p = 1.0
    for i in 1:nq
        qinds = qinds_func(i, nq)
        p *= grouped_marginal(running_psum, i, x[i], qinds)
        running_psum = project(running_psum, i, x[i]) # TODO: Can I use the trace_project
    end
    return p
end

function grouped_prob(psum, x::Integer; qinds_func=qinds_base_func)
    nq = psum.nqubits
    return grouped_prob(psum, BitVector(digits(x, base=2, pad=nq)); qinds_func = qinds_func)
end

function bayes_unnormalized_marginal(psum::PauliSum{TT, CT}, qind::Integer, bit::Bool) where {TT, CT}
    a = getcoeff(psum, :I, qind)
    b = getcoeff(psum, :Z, qind)

    return abs(a + (1 - 2 * bit) * b)
end

function bayes_marginal(psum::PauliSum{TT, CT}, qind::Integer, bit::Bool)
    p0 = bayes_unnormalized_marginal(psum, qind, false)
    p1 = bayes_unnormalized_marginal(psum, qind, true)

    return (bit == 0 ? p0 : p1) / (p0 + p1)
end

function bayes_prob(psum::PauliSum{TT, CT}, x::BitVector) where {TT, CT}
    nq = psum.nqubits
    running_psum = psum
    p = 1.0
    for i in 1:nq
        p *= bayes_marginal(running_psum, i, x[i])
        running_psum = trace_project(running_psum, i, x[i])
    end
    return p
end

function bayes_prob(psum, x::Integer)
    return prob(psum, BitVector(digits(x, base=2, pad=nq)))
end
