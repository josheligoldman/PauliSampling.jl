using Test
using PauliSampling

@testset "Normalizing factor consistency" begin
    num_set_bits = 3
    set_bits = rand(Bool, num_set_bits) |> BitVector
    nq = 6
    eps = 1e-6

    for _ in 1:10
        psum = PauliSampling.random_psum(nq)
        for i in 1:num_set_bits
            psum = PauliSampling.apply_projector(psum, i, set_bits[i])
        end

        x1 = PauliSampling.normalizing_factor(psum)
        x2 = PauliSampling.faster_normalizing_factor(psum)
        x3 = coset_normalizing_factor(psum; num_set_bits=num_set_bits, set_bits=set_bits)

        @test isapprox(x1, x2; atol=eps)
        @test isapprox(x2, x3; atol=eps)
        @test isapprox(x3, x1; atol=eps)
    end
end