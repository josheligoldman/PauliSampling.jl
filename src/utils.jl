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

function qinds_in_order(i, num, nq)
    return i:min(i + (num-1), nq)
end

function qinds_random(i, num, nq)
    idx_range = i:nq
    return k ≤ length(idx_range) ? shuffle(idx_range)[1:num] : collect(idx_range)
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

function get_dist(psum, prob_fn)
    nq = psum.nqubits
    dist = Vector{Float64}(undef, 2^nq)
    for i in 0:(2^nq - 1)
        dist[i+1] = prob_fn(psum, i)
    end
    return dist
end

function tvd(p, q)
    return 0.5 * sum(abs.(p .- q))
end

function kl_div(p, q)
    res = 0.0
    for (i, pi) in enumerate(p)
        qi = q[i]
        res += pi * log(pi/qi)
    end
    return res
end

@inline function get_bit(a::Unsigned, i::Integer)::Bool
    return ((a >> (i - 1)) & 0x01) == 1
end

@inline function set_bit(a::T, i::Integer, b::Bool)::T where {T <: Unsigned}
    mask = one(T) << (i-1)
    if b
        return a | mask
    else
        return a & ~mask
    end
end