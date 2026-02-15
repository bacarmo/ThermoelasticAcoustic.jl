using TestItemRunner

@run_package_tests verbose = true

@testitem "my_sum: type stability" begin
    using ThermoelasticAcoustic: my_sum

    integers = [1, 2, 3]
    floats = [1.0, 2.0, 3.0]
    float32s = Float32[1, 2, 3]

    @inferred my_sum(integers)
    @inferred my_sum(floats)
    @inferred my_sum(float32s)

    @test my_sum(integers) isa Int
    @test my_sum(floats) isa Float64
    @test my_sum(float32s) isa Float32
end

@testitem "my_sum: correctness" begin
    using ThermoelasticAcoustic: my_sum

    @test my_sum(1:10) == 55
    @test my_sum([1, 2, 3, 4]) == 10
    @test my_sum([1.0, 2.0, 3.0]) ≈ 6.0
end

@testitem "my_sum: edge cases" begin
    using ThermoelasticAcoustic: my_sum

    @test my_sum(Int[]) == 0
    @test my_sum(1:0) == 0
    @test my_sum(Float64[]) == 0.0
end