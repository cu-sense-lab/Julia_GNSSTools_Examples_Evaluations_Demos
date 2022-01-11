using Statistics
using LinearAlgebra
using Optim
using QuadGK


function init_H(k₁::Number, k₂::Number)
    H(s) = (k₁*s + k₂) / (s^2 + k₁*s + k₂)
    return H
end


function init_H(k₁::Number, k₂::Number, k₃::Number)
    H(s) = (k₁*s^2 + k₂*s + k₃) / (s^3 + k₁*s^2 + k₂*s + k₃)
    return H
end


function init_H2(T, k₁::Number, k₂::Number, k₃::Number)
    α = k₁
    β = k₂
    γ = k₃
    numerator_(z) = (7*T^2*γ + 9*T*β + 6*α)*z^2 +
                    (-2*T^2*γ - 12*T*β - 12*α)*z +
                    T^2*γ + 3*T*β + 6*α
    denominator_(z) = 6*z^3 +
                      (7*T^2*γ + 9*T*β + 6*α - 18)*z^2 +
                      (-2T^2*γ - 12*T*β - 12*α + 18)*z +
                      T^2*γ + 3*T*β + 6*α - 6
    H(z) = numerator_(z) / denominator_(z)
end


# function get_B(H, f=1e30)
#     h(f) = abs2(H(2π*f*1im))
#     return real(quadgk(h, 0, f)[1]) / abs2(H(0))
# end


function get_B(H, f=1e6)
    s = big(2π*f*1im)
    h(s) = H(s) * H(-s) / (2π*1im)
    return real(quadgk(h, -s, s)[1]) / abs2(H(0)) / 2
end


function get_B_from_phase_noise(σᵩ², CN0, T, h₀=0, h₋₁=0, h₋₂=0)
    f(x) = abs(phase_noise_variance(x[1], CN0, T, h₀, h₋₁, h₋₂) - σᵩ²)
    initial_B = BigFloat.([15])
    B = optimize(f, initial_B, LBFGS()).minimizer[1]
    return B
end


function calc_ωₙ(B)
    return B/0.7845
end


function calc_B_from_ωₙ(ωₙ, aₙ, bₙ)
    B = ωₙ*(aₙ*bₙ^2 + aₙ^2 - bₙ) / (aₙ*bₙ - 1) / 4
    return B
end


function calc_L_2state(T, B, damping=0.707)
    ωₙ = B/0.53
    α = 2*damping*ωₙ*T - 3*ωₙ^2*T^2/2
    β = ωₙ^2*T
    return [α, β]
end


function calc_L(T, ωₙ, aₙ=1.1, bₙ=2.4)
    α = (11*ωₙ^3*T^3 - 9*aₙ*ωₙ^2*T^2 + 6*bₙ*ωₙ*T)/6
    β = -2*ωₙ^3*T^2 + aₙ*ωₙ^2*T
    γ = ωₙ^3*T
    return [α, β, γ]
end


function get_B_from_rong(T, k₁, k₂, k₃; initial_aₙ=1.1, initial_bₙ=2.4,
                         initial_B=15)
    function f(x)
        α, β, γ = calc_L(T, x[1], x[2], x[3])
        δα = abs2(k₁ - α)
        δβ = abs2(k₂ - β)
        δγ = abs2(k₃ - γ)
        err = δα + δβ + δγ
        return err
    end
    initial_x = BigFloat.([calc_ωₙ(initial_B), initial_aₙ, initial_bₙ])
    ωₙ, aₙ, bₙ = optimize(f, initial_x, LBFGS()).minimizer
    B = calc_B_from_ωₙ(ωₙ, aₙ, bₙ)
    return Float64(B)
end


function calc_B(k₁, k₂)
    return (k₁ + k₂/k₁)/4
end


function calc_B(k₁, k₂, k₃)
    return (k₁ + k₂^2/(k₂*k₁ - k₃))/4
end


function calc_B(T, k₁, k₂, k₃)
    m₀ = 2π*k₃/T
    m₁ = 2π*k₃ + 2π*k₂/T
    m₂ = T*π*k₃/3 + π*k₂ + k₁*T
    B = m₂/4 + m₁^2 / (m₁*m₂ - m₀) / 4
    return B
end
