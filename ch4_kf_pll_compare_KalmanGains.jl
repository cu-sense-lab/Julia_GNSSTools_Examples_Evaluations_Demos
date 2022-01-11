using GNSSTools
using MatrixEquations
using LinearAlgebra
using PyPlot
pygui(true)
include("kalman_filter_gain_error_functions.jl")

h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
# h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
h₀ = h_parms[3]
h₋₁ = h_parms[2]
h₋₂ = h_parms[1]

T = 1e-3
qₐ = 1
f_L = L1_freq
state_num = 3
CN0 = 45
B = 15

A = calcA(T, state_num)
C = calcC(T, state_num)
# R = [phase_noise_variance(20, T)]
R = [55*π/180] .^2
Q = calcQ(T, h₀, h₋₂, qₐ, f_L, state_num)
P, CLSEIG, F = ared(A', C', R[1], Q)
K = (P*C')*pinv(C*P*C' + R)
pᵩ = sqrt((C*P*C')[1]) * 180/π

σ²ᵩ = phase_noise_variance(CN0, T)
qₐ = 0.2^2
ω³₀ = 2π*f_L*sqrt(qₐ/(σ²ᵩ*T))/GNSSTools.c
ω₀ = (ω³₀)^(1/3)
ω²₀ = ω₀^2
K2 = [2*ω₀, 2*ω²₀, ω³₀]
B2 = calc_B(K2...)

# h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
h₀ = h_parms[3]
h₋₁ = h_parms[2]
h₋₂ = h_parms[1]
qₐ = 0.2^2
CN0s = [25, 35, 45]
Ts = Array(range(1.0, 20.0, step=1.0))*1e-3
fig = figure(figsize=(7.5, 4))
ax1 = fig.add_subplot(1,2,1)
xlabel("Integration Time (ms)")
ylabel("Phase Uncertainty (deg)")
ax2 = fig.add_subplot(1,2,2)
xlabel("Integration Time (ms)")
ylabel("Equivalent Bandwidth (Hz)")
for CN0 in CN0s
    σᵩ = zeros(length(Ts))
    Bs = zeros(length(Ts))
    for i in 1:length(Ts)
        T = Ts[i]
        A = calcA(T, state_num)
        C = calcC(T, state_num)
        R = [phase_noise_variance(CN0, T)]
        Q = calcQ(T, h₀, h₋₂, qₐ, f_L, state_num)
        P, CLSEIG, F = ared(A', C', R[1], Q)
        # pᵩ = sqrt((C*P*C')[1])
        pᵩ = sqrt(P[1])
        K = (P*C')*pinv(C*P*C' + R)
        σᵩ[i] = pᵩ*180/π
        Bs[i] = calc_B(K...)
        # Bs[i] = get_B(init_H(K...))
    end
    ax1.plot(Ts.*1000, σᵩ, label="$(CN0) dB⋅Hz")
    ax2.plot(Ts.*1000, Bs, label="$(CN0) dB⋅Hz")
end
legend()
subplots_adjust(wspace=0.3, bottom=0.15, left=0.08, right=0.92, top=0.92)


# h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
h₀ = h_parms[3]
h₋₁ = h_parms[2]
h₋₂ = h_parms[1]
qₐ = 0.05^2
CN0s = [10, 15, 20, 25]
Ts = Array(range(20, 140, step=1.0))*1e-3
fig = figure()
for CN0 in CN0s
    σᵩ = zeros(length(Ts))
    for i in 1:length(Ts)
        T = Ts[i]
        A = calcA(T, state_num)
        C = calcC(T, state_num)
        R = [phase_noise_variance(CN0, T)]
        Q = calcQ(T, h₀, h₋₂, qₐ, f_L, state_num)
        P, CLSEIG, F = ared(A', C', R[1], Q)
        K = (P*C')*pinv(C*P*C' + R)
        σᵩ[i] = sqrt(P[1])*180/π
    end
    plot(Ts.*1000, σᵩ, label="$(CN0) dB⋅Hz")
end
xlabel("Integration Time (ms)")
ylabel("Phase Uncertainty (deg)")
legend()
