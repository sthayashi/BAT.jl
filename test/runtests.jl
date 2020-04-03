# This file is a part of BAT.jl, licensed under the MIT License (MIT).

using Test

Test.@testset "Package BAT" begin
    include("integration/test_bat_integrate.jl")
end
