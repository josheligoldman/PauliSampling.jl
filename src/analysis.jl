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

function shannon_entropy(prob)
    entropy = 0.0
    for p in prob
        if p > 0  # avoid log(0)
            entropy -= p * log2(p)
        end
    end
    return entropy
end



function avg_metrics_for_weight(nq, max_weight, M)
    tvd_coset = 0.0
    tvd_bayes = 0.0
    kl_coset  = 0.0
    kl_bayes  = 0.0

    circuit_nq  = build_circuit(nq)
    nparams  = countparameters(circuit_nq)
    init_psum = zero_state(nq)

    for _ in 1:M
        thetas = randn(nparams) * 0.5
        exact_psum = propagate(circuit_nq, init_psum, thetas)
        trunc_psum = propagate(circuit_nq, init_psum, thetas; max_weight=max_weight)

        exact_dist = PauliSampling.get_dist(exact_psum,  PauliSampling.approximate_prob)
        coset_dist = PauliSampling.get_dist(trunc_psum, PauliSampling.projection_prob)
        bayes_dist = PauliSampling.get_dist(trunc_psum, PauliSampling.approximate_prob)

        tvd_coset += PauliSampling.tvd(exact_dist, coset_dist)
        kl_coset += PauliSampling.kl_div(exact_dist, coset_dist)
        tvd_bayes += PauliSampling.tvd(exact_dist, bayes_dist)
        kl_bayes += PauliSampling.kl_div(exact_dist, bayes_dist)
    end
    return tvd_coset/M, kl_coset/M, tvd_bayes/M, kl_bayes/M
end

function boolean_fourier_dist_diff(nq, p, q)
    dist_diff = p - q
    fourier_coeffs = naive_fwht(dist_diff)
    bins = zeros(nq + 1)
    for i in 1:2^nq
        bins[count_ones(i-1) + 1] += abs(fourier_coeffs[i])
    end

    return bins
end

function naive_fwht!(x::Vector{Float64})
    n = length(x)
    @assert ispow2(n) "Length must be a power of 2"
    h = 1
    while h < n
        for i in 1:2h:n
            for j in i:i+h-1
                u = x[j]
                v = x[j + h]
                x[j] = u + v
                x[j + h] = u - v
            end
        end
        h *= 2
    end
    return x
end

function naive_fwht(x::Vector{Float64})
    return naive_fwht!(copy(x))
end

function compare_dists_boolean(nq; max_weight=4)
    init_psum = zero_state(nq)
    circuit  = build_circuit(nq)
    nparams  = countparameters(circuit)
    thetas = randn(nparams) * 0.5
    exact_psum = propagate(circuit, init_psum, thetas)
    trunc_psum = propagate(circuit, init_psum, thetas; max_weight=max_weight)
    
    exact_dist = get_dist(exact_psum, approximate_prob)
    coset_dist = get_dist(trunc_psum, projection_prob)
    bayes_dist = get_dist(trunc_psum, approximate_prob)

    exact_coset_bins = boolean_fourier_dist_diff(nq, exact_dist, coset_dist)
    exact_bayes_bins = boolean_fourier_dist_diff(nq, exact_dist, bayes_dist)

    x_bins = 0:nq
    bin_plt = plot(x_bins, exact_coset_bins; label="Coset vs Exact", lw=2, marker=:circle)
    plot!(x_bins, exact_bayes_bins; label="Bayes vs Exact", lw=2, marker=:square)
end



function bit_marginals(nq, dist)
    marg = zeros(nq)
    for i in 1:length(dist)
        bits = digits(i - 1, base=2, pad=nq)
        for j in 1:nq
            marg[j] += bits[j] * dist[i]
        end
    end
    return marg
end