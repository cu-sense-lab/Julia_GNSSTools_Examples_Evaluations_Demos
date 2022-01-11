using GNSSTools
using ProgressMeter
using Distributions
using JLD
using PyPlot
pygui(true)
include("kalman_filter_gain_error_functions.jl")

tcxo_h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
ocxo_h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
h_parms = tcxo_h_parms
h₋₂ = h_parms[1]
h₋₁ = h_parms[2]
h₀ = h_parms[3]
sig_freq = L1_freq
t_length = 1  # s
phi_init = 0  # rad
f_d_init = 1000  # Hz
fd_rate_init = -0.7  # Hz/s


#----------------------------------SIMULATION----------------------------------#

σω = 1
q_a = 0.1
T = 1e-3
M = floor(Int, t_length/T)
t = calctvector(M, 1/T)

A = calcA(T, 3)
C = calcC(T, 3)
Q = calcQ(T, h₀, h₋₂, q_a, sig_freq, 3)
R = [phase_noise_variance(CN0, T)]

x⁺_init = [phi_init; 2π*f_d_init; 2π*fd_rate_init]
x⁺ = deepcopy(x⁺_init)
x_truth_noiseless = zeros(3, M)
ϕ_truth_noiseless = zeros(M)
x_truth = zeros(3, M)
ϕ_truth = zeros(M)
Δϕ_truth = zeros(M)

# Generate truth
x_truth_noiseless[:,1] = x⁺_init
x⁺ = deepcopy(x⁺_init)
for i in 2:M
    x⁺ = A*x⁺
    x_truth_noiseless[:,i] = x⁺
    ϕ_truth_noiseless[i] = (C*x⁺)[1]
end

# Simulate states
x_truth[:,1] = x⁺_init
x⁺ = deepcopy(x⁺_init)
for i in 2:M
    x⁺ = A*x⁺ + rand(MvNormal(Q))
    x_truth[:,i] = x⁺
    ϕ_truth[i] = (C*x⁺ + rand(MvNormal(R)))[1]
    if i > 1
        Δϕ_truth[i] = (C*(x_truth_noiseless[:,i] .- x_truth[:,i]))[1]
    end
end

fig = figure(figsize=(8, 5))

ax1 = fig.add_subplot(3,2,1)
ax1.plot(t, x_truth[1,:], "k-")
xlabel("Time (s)")
ylabel("ϕ State (rad)")

ax2 = fig.add_subplot(3,2,2)
ax2.plot(t, ϕ_truth, "k-")
xlabel("Time (s)")
ylabel("ϕ Meas. (rad)")

ax3 = fig.add_subplot(3,2,3)
ax3.plot(t, x_truth[2,:]./2π, "k-")
xlabel("Time (s)")
ylabel("f_d (Hz)")

ax4 = fig.add_subplot(3,2,4)
ax4.plot(t, Δϕ_truth.*(180/π), "k-")
xlabel("Time (s)")
ylabel("Δϕ Meas. (degrees)")

ax5 = fig.add_subplot(3,2,5)
ax5.plot(t, x_truth[3,:]./2π, "k-")
xlabel("Time (s)")
ylabel("f_d rate (Hz/s)")

ax6 = fig.add_subplot(3,2,6)
ax6.plot(t, (ϕ_truth_noiseless .- ϕ_truth).*(180/π), "k-")
# ax6.plot(t, ϕ_truth_noiseless.*(180/π), "k:")
xlabel("Time (s)")
ylabel("ϕ Truth (rad)")

subplots_adjust(wspace=0.4, hspace=0.6, bottom=0.15, left=0.12, right=0.88, 
                top=0.93)

#------------------------------------FILTER------------------------------------#

CN0 = 45
σω = 1
q_a = 0.1
T = 1e-3
state_num = 3

A = calcA(T, state_num)
C = calcC(T, state_num)
Q = calcQ(T, h₀, h₋₂, q_a, sig_freq, state_num)
R = [phase_noise_variance(CN0, T)]

if state_num == 3
    P⁺ = [phase_noise_variance(B_init, CN0, T, h₀, h₋₁, h₋₂) 0 0;
          0 (2π/T)^2 0;
          0 0 σω^2]
    x_init = [phi_init; 2π*f_d_init; 2π*fd_rate_init]# + rand(MvNormal(P⁺))
else
    P⁺ = [phase_noise_variance(B_init, CN0, T, h₀, h₋₁, h₋₂) 0;
          0 (2π/T)^2]
    x_init = [phi_init; 2π*f_d_init]# + rand(MvNormal(P⁺))
end

K = zeros(state_num, M)
P = zeros(state_num, M)
x = zeros(state_num, M)
x[:,1] = x_init
P[:,1] = diag(P⁺)
pᵩ = zeros(M)
x⁺ = deepcopy(x_init)
for i in 2:M
    x⁻ = A*x⁺
    P⁻ = A*P⁺*A' + Q
    Δϕ = Δϕ_truth[i]
    Kᵢ = (P⁻*C')*pinv(C*P⁻*C' + R)
    P⁺ = (I - Kᵢ*C)*P⁻
    x⁺ = x⁻ + Kᵢ*[Δϕ]
    K[:,i] = Kᵢ
    P[:,i] = diag(P⁺)
    x[:,i] = x⁺
    pᵩ[i] = (C*P⁺*C')[1]
end

fig = figure(figsize=(8, 5))

ax1 = fig.add_subplot(3,2,1)
ax1.plot(t, x_truth_noiseless[1,:] .- x[1,:], "k-")
# ax1.plot(t, x_truth_noiseless[1,:], "k:")
xlabel("Time (s)")
ylabel("ϕ State (rad)")

ax2 = fig.add_subplot(3,2,2)
ax2.plot(t, sqrt.(pᵩ).*(180/π), "k-")
xlabel("Time (s)")
ylabel("pᵩ (degrees)")

ax3 = fig.add_subplot(3,2,3)
ax3.plot(t, x_truth_noiseless[2,:]./2π .- x[2,:]./2π, "k-")
# ax3.plot(t, x_truth_noiseless[2,:]./2π, "k:")
xlabel("Time (s)")
ylabel("f_d (Hz)")

ax4 = fig.add_subplot(3,2,4)
ax4.plot(t, Δϕ_truth.*(180/π), "k-")
xlabel("Time (s)")
ylabel("Δϕ Meas. (degrees)")

if state_num == 3
    ax5 = fig.add_subplot(3,2,5)
    ax5.plot(t, x_truth_noiseless[3,:]./2π .- x[3,:]./2π, "k-")
    # ax1.plot(t, x_truth_noiseless[3,:]./2π, "k:")
    xlabel("Time (s)")
    ylabel("f_d rate (Hz/s)")
end

ax6 = fig.add_subplot(3,2,6)
ax6.plot(t, K[1,:], label="k₁")
ax6.plot(t, K[2,:], label="k₂")
if state_num == 3
    ax6.plot(t, K[3,:], label="k₃")
end
xlabel("Time (s)")
ylabel("KF Gain")

subplots_adjust(wspace=0.4, hspace=0.6, bottom=0.15, left=0.12, right=0.88, 
                top=0.93)



##############

function test_KF(CN0, T)
    f_s = 2.5e6
    f_if = 0
    sig_freq = L1_freq
    B = 2*l1ca_chipping_rate
    Tsys = 535
    nADC = 8
    t_length = 1
    tcxo_h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
    h₋₂ = tcxo_h_parms[1]
    h₋₁ = tcxo_h_parms[2]
    h₀ = tcxo_h_parms[3]
    f_d0 = 1000
    f_dd0 = 0
    phi0 = 1.2
    state_num = 3
    qₐ = 1

    carrier_amp = sqrt(2*GNSSTools.k*Tsys)*10^(CN0/20)  # for sinusoid
    noise_amp = sqrt(GNSSTools.k*B*Tsys)
    adc_scale = 2^(nADC-1)
    N = floor(Int, t_length*f_s)
    t = calctvector(N, f_s)
    phase_noise = real.(generate_phase_noise(t_length, f_s, sig_freq, tcxo_h_parms))
    thermal_noise = noise_amp .* randn(ComplexF64, N)
    signal = carrier_amp.*cis.((2π).*(f_if .+ f_d0 .+ 0.5 .* f_dd0.*t).*t .+ phi0)
    phi_true = (2π).*(f_if .+ f_d0 .+ 0.5 .* f_dd0.*t.^2).*t .+ phi0
    data = signal.*cis.(phase_noise) .+ thermal_noise
    # data = signal
    sigmax = sqrt(maximum(abs2, data))
    # data = round.(data.*(adc_scale/sigmax))
    for i in 1:N
        data[i] = round(data[i]*adc_scale/sigmax)
        # Check that the max positive value for I/Q channels is 
        # (2^(nADC - 1) - 1)
        if real(data[i]) == adc_scale
            data[i] = data[i] - 1 
        end
        if imag(data[i]) == adc_scale
            data[i] = data[i] - 1im
        end
    end

    M = floor(Int, t_length/T)
    M_over_2 = floor(Int, M/2)
    Mₙ = floor(Int, T*f_s)
    σ²ᵩ_unfiltered = phase_noise_variance(CN0, T)
    A = calcA(T, state_num)
    C = calcC(T, state_num)
    R = [σ²ᵩ_unfiltered]
    Q = calcQ(T, h₀, h₋₂, qₐ, sig_freq, state_num)
    x = zeros(state_num, M)
    Δx = zeros(state_num, M)
    P = zeros(state_num, M)
    K = zeros(state_num, M)
    Δϕs = zeros(M)
    ZP = zeros(ComplexF64, M)
    integrated_phase_noise = zeros(M)
    if state_num == 3
        P⁻ = diagm([sqrt(σ²ᵩ_unfiltered), 2π/(2T), 2π*5]).^2
        x⁻ = [phi0, 2π*f_d0, 2π*f_dd0] .+ rand(MvNormal(P⁻))
    elseif state_num == 2
        P⁻ = diagm([sqrt(σ²ᵩ_unfiltered), 2π/(2T)]).^2
        x⁻ = [phi0, 2π*f_d0] .+ rand(MvNormal(P⁻))
    else
        @error("Invalid state_num value of $(state_num). Only 2 and 3 are valid.")
    end
    x⁻ = x⁻ .+ rand(MvNormal(P⁻))
    track_t = calctvector(M, 1/T)

    for i in 1:M
        if state_num == 3
            ϕ, ω, ω_dot = x⁻
        else
            ϕ, ω = x⁻
            ω_dot = 0
        end
        subdata = data[(i-1)*Mₙ+1:i*Mₙ]
        replica = cis.(-(2π.*f_if.*t[1:Mₙ] .+ ω.*t[1:Mₙ] .+ ω_dot.*t[1:Mₙ].^2 ./ 2 .+ ϕ))
        zp = sum(subdata.*replica)/Mₙ
        # Δϕ = atan(imag(zp)/real(zp))
        Δϕ = atan(imag(zp), real(zp))
        Kᵢ = (P⁻*C')*pinv(C*P⁻*C' + R)
        P⁺ = (I - Kᵢ*C)*P⁻*(I - Kᵢ*C)' + Kᵢ*R*Kᵢ'
        correction = Kᵢ.*Δϕ
        x⁺ = x⁻ + correction
        x[:,i] = x⁺
        Δx[:,i] = correction
        P[:,i] = diag(P⁺)
        K[:,i] = Kᵢ
        Δϕs[i] = Δϕ
        ZP[i] = zp
        sub_phase_noise = phase_noise[(i-1)*Mₙ+1:i*Mₙ]
        integrated_phase_noise[i] = mean(sub_phase_noise)
        x⁻ = A*x⁺
    end
    return (track_t, x, P, K, Δϕs, Δx, ZP, integrated_phase_noise, 
            N, M, Mₙ, M_over_2, σ²ᵩ_unfiltered, phi_true)
end

CN0 = 45
T = 5e-3

track_t, x, P, K, Δϕs, Δx, ZP, integrated_phase_noise, 
N, M, Mₙ, M_over_2, σ²ᵩ_unfiltered, phi_true = test_KF(CN0, T)

fig = figure()
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(track_t, Δϕs[:].*(180/π), "k.", label=L"\Delta\phi")
ax1.plot(track_t, Δx[1,:].*(180/π), label=L"\Delta\hat{\phi}")
ax1.plot(track_t, sqrt.(P[1,:]).*(180/π), 
         linestyle=":", color="r")
ax1.plot(track_t, -sqrt.(P[1,:]).*(180/π), linestyle=":", color="r")
unfiltered_sigma_phi = std(Δϕs[M_over_2:end])
filtered_sigma_phi = std(Δx[1,M_over_2:end])
ax1.plot([track_t[1], track_t[end]], [unfiltered_sigma_phi, unfiltered_sigma_phi].*(180/π), 
          linestyle=":", color="g")
ax1.plot([track_t[1], track_t[end]], 
         [-unfiltered_sigma_phi, -unfiltered_sigma_phi].*(180/π),
         linestyle=":", color="g")
ax1.plot([track_t[1], track_t[end]], 
         [filtered_sigma_phi, filtered_sigma_phi].*(180/π), 
          linestyle=":", color="gray")
ax1.plot([track_t[1], track_t[end]], 
         [-filtered_sigma_phi, -filtered_sigma_phi].*(180/π),
         linestyle=":", color="gray")
ax1.plot([track_t[1], track_t[end]], 
         [sqrt(σ²ᵩ_unfiltered), sqrt(σ²ᵩ_unfiltered)].*(180/π), 
          linestyle=":", color="orange")
ax1.plot([track_t[1], track_t[end]], 
         [-sqrt(σ²ᵩ_unfiltered), -sqrt(σ²ᵩ_unfiltered)].*(180/π),
         linestyle=":", color="orange")
legend(ncol=1)
xlabel("Time (seconds)")
ylabel("Phase (deg)")

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(track_t, (phi_true[1:Mₙ:N] .+ integrated_phase_noise .- x[1,:]).*(180/π), 
         label=L"\phi + \phi_{oscillator} - \hat{\phi}")
ax2.plot(track_t, sqrt.(P[1,:]).*(180/π), label=L"\sqrt{P_{1,1}}", 
         linestyle=":", color="r")
ax2.plot(track_t, -sqrt.(P[1,:]).*(180/π), linestyle=":", color="r")
unfiltered_sigma_phi = std(Δϕs[M_over_2:end])
filtered_sigma_phi = std(Δx[1,M_over_2:end])
ax2.plot([track_t[1], track_t[end]], [unfiltered_sigma_phi, unfiltered_sigma_phi].*(180/π), 
          label=L"\sigma_{\Delta\phi}", linestyle=":", color="g")
ax2.plot([track_t[1], track_t[end]], 
         [-unfiltered_sigma_phi, -unfiltered_sigma_phi].*(180/π),
         linestyle=":", color="g")
ax2.plot([track_t[1], track_t[end]], 
         [filtered_sigma_phi, filtered_sigma_phi].*(180/π), 
          label=L"\sigma_{\Delta\hat{\phi}}", linestyle=":", color="gray")
ax2.plot([track_t[1], track_t[end]], 
         [-filtered_sigma_phi, -filtered_sigma_phi].*(180/π),
         linestyle=":", color="gray")
ax2.plot([track_t[1], track_t[end]], 
         [sqrt(σ²ᵩ_unfiltered), sqrt(σ²ᵩ_unfiltered)].*(180/π), 
          label=L"\sqrt{\sigma^2_\phi}}", linestyle=":", color="orange")
ax2.plot([track_t[1], track_t[end]], 
         [-sqrt(σ²ᵩ_unfiltered), -sqrt(σ²ᵩ_unfiltered)].*(180/π),
         linestyle=":", color="orange")
legend()
xlabel("Time (seconds)")
ylabel("Phase (deg)")

ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(track_t, real.(ZP), label="real(ZP)")
ax3.plot(track_t, imag.(ZP), label="imag(ZP)")
legend()
xlabel("Time (seconds)")
ylabel("ZP")


function test_KF_errs(CN0, T, iterations)
    P₁₁ = 0.
    σ²_Δϕ_unfilt = 0.
    σ²_Δϕ_filt = 0.
    σ²ᵩ_unfiltered = phase_noise_variance(CN0, T)
    for iteration in 1:iterations
        track_t, x, P, K, Δϕs, Δx, ZP, integrated_phase_noise, 
        N, M, Mₙ, M_over_2, σ²ᵩ_unfiltered, phi_true = test_KF(CN0, T)
        σ²_Δϕ_unfilt += var(Δϕs[M_over_2:end])
        σ²_Δϕ_filt += var(Δx[1,M_over_2:end])
        P₁₁ += mean(P[1,M_over_2:end])
    end
    P₁₁_sqrt = sqrt(P₁₁/iterations)
    σ_Δϕ_unfilt = sqrt(σ²_Δϕ_unfilt/iterations)
    σ_Δϕ_filt = sqrt(σ²_Δϕ_filt/iterations)
    σᵩ_unfiltered = sqrt(σ²ᵩ_unfiltered)
    return (P₁₁_sqrt, σ_Δϕ_unfilt, σ_Δϕ_filt, σᵩ_unfiltered)
end

CN0s = Array(range(20, 50, step=1))
Ts = [1e-3, 2e-3]
# Ts = [1e-3]
iterations = 5

P₁₁_sqrts = []
σ_Δϕ_unfilts = []
σ_Δϕ_filts = []
σᵩ_unfiltereds = []

for j in 1:length(Ts)
    T = Ts[j]
    N = length(CN0s)
    P₁₁_sqrts_ = zeros(N)
    σ_Δϕ_unfilts_ = zeros(N)
    σ_Δϕ_filts_ = zeros(N)
    σᵩ_unfiltereds_ = zeros(N)
    p = Progress(N, 1, "Processing...")
    for i in 1:N 
        CN0 = CN0s[i]
        P₁₁_sqrt, σ_Δϕ_unfilt, σ_Δϕ_filt, σᵩ_unfiltered = test_KF_errs(CN0, T, iterations)
        # println([P₁₁_sqrt, σ_Δϕ_unfilt, σ_Δϕ_filt, σᵩ_unfiltered])
        P₁₁_sqrts_[i] = P₁₁_sqrt
        σ_Δϕ_unfilts_[i] = σ_Δϕ_unfilt
        σ_Δϕ_filts_[i] = σ_Δϕ_filt
        σᵩ_unfiltereds_[i] = σᵩ_unfiltered
        next!(p)
    end
    push!(P₁₁_sqrts, P₁₁_sqrts_)
    push!(σ_Δϕ_unfilts, σ_Δϕ_unfilts_)
    push!(σ_Δϕ_filts, σ_Δϕ_filts_)
    push!(σᵩ_unfiltereds, σᵩ_unfiltereds_)
end


fig = figure(figsize=(7.5, 5))
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(CN0s, P₁₁_sqrts[1].*(180/π), label=L"\sqrt{P_{11}}")
ax1.plot(CN0s, σ_Δϕ_filts[1].*(180/π), label=L"\sigma_{\Delta\hat{\phi}}")
hlines(y=10log10(15), xmin=CN0s[1], xmax=CN0s[end], color="grey", 
       linestyle=":", label=string(L"15^\circ", " Tracking Limit"))
ax1.set_yscale("log")
legend()
xlabel("C/N₀ (dB⋅Hz)")
ylabel("Carrier Phase Error (deg)")

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(CN0s, σᵩ_unfiltereds[1].*(180/π), label=L"\sigma_{\phi}")
ax2.plot(CN0s, σ_Δϕ_unfilts[1].*(180/π), label=L"\sigma_{\Delta\phi}")
ax2.set_yscale("log")
legend()
xlabel("C/N₀ (dB⋅Hz)")
ylabel("Carrier Phase Error (deg)")


ax3 = fig.add_subplot(2, 2, 3)
ax3.plot(CN0s, P₁₁_sqrts[2].*(180/π), label=L"\sqrt{P_{11}}")
ax3.plot(CN0s, σ_Δϕ_filts[2].*(180/π), label=L"\sigma_{\Delta\hat{\phi}}")
hlines(y=10log10(15), xmin=CN0s[1], xmax=CN0s[end], color="grey", 
       linestyle=":", label=string(L"15^\circ", " Tracking Limit"))
ax3.set_yscale("log")
legend()
xlabel("C/N₀ (dB⋅Hz)")
ylabel("Carrier Phase Error (deg)")

ax4 = fig.add_subplot(2, 2, 4)
ax4.plot(CN0s, σᵩ_unfiltereds[2].*(180/π), label=L"\sigma_{\phi}")
ax4.plot(CN0s, σ_Δϕ_unfilts[2].*(180/π), label=L"\sigma_{\Delta\phi}")
ax4.set_yscale("log")
legend()
xlabel("C/N₀ (dB⋅Hz)")
ylabel("Carrier Phase Error (deg)")

subplots_adjust(wspace=0.3, bottom=0.15, left=0.08, right=0.92, top=0.92)

