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