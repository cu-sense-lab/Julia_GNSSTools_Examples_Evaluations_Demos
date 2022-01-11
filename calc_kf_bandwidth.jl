using QuadGK
using Statistics
using GNSSTools
using ProgressMeter
using JLD
using Random
using Optim
using PyPlot
pygui(true)
include("kalman_filter_gain_error_functions.jl")


# k₁ = 0.017299
# k₂ = 0.024087
# B2_truth = (k₁ + k₂/k₁)/4
# H2 = init_H(k₁, k₂)
# B2 = get_B(H2)

# H3 = init_H(0.5, 2.0, 4)
# B3 = get_B(H3)


# h₀ = 1e-19
# h₋₂ = 2e-20
# T = 0.001  # s
# qₐ = 0.1
# CN0 = 30 # dB⋅Hz
# sig_freq = L1_freq
# state_num = 2
# # R = [1e-4]
# R = [1/(2*10^(CN0/10)*T)*(1 + 1/(2*10^(CN0/10)*T))]
# A = calcA(T, state_num)
# C = calcC(T, state_num)
# Q = calcQ(T, h₀, h₋₂, qₐ, sig_freq, state_num)
# K = dkalman(A, C, Q, GNSSTools.Diagonal(R))
# if length(K) == 2
#     k1 = K[1]
#     k2 = K[2]
#     B_truth = (k1 + k2/k1)/4
#     B = get_B(init_H(k1, k2))
#     println([B, B_truth, B_truth-B])
# else
#     k1 = K[1]
#     k2 = K[2]
#     k3 = K[3]
#     B = get_B(init_H(k1, k2, k3))
#     println(B)
# end

# ############################

# t_length = 1
# f_s = 5e6
# # f_s = 2.048e6
# # f_s = 25e6
# # h_parms = h_parms_tcxo[1]
# h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
# # h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO

# sigtype = define_l1ca_code_type(; sig_freq=1*L1_freq)
# # sigtype = define_l5_code_type(t_length; sig_freq=1*L5_freq, channel="both")
# # sigtype = definesignaltype(definecodetype(l5i_codes, L5_chipping_rate), 
#                         #    L5_freq, "I")
# # signal = definesignal(sigtype, 5e6, t_length; 
# #                       prn=1, code_start_idx=1 000, f_d=163000, fd_rate=449, 
# #                       CN0=45, include_phase_noise=false, 
# #                       include_thermal_noise=true, receiver_h_parms=h_parms);
# signal = definesignal(sigtype, f_s, t_length; 
#                       prn=1, code_start_idx=1000, f_d=1000, fd_rate=0, f_if=1.25e6,
#                       CN0=25, include_phase_noise=true, 
#                       include_thermal_noise=true, include_adc=true, 
#                       include_carrier=true, receiver_h_parms=h_parms)
# generatesignal!(signal);
# # signal.noexp = false
# # generatereplica!(signal);

# acqresults, trackresults, corr_result, SNR_est = process(signal, 
#                                                          sigtype, 
#                                                          1, "I"; fd_range=5000, 
#                                                          fine_acq_method=:fft, 
#                                                          return_corrresult=true, 
#                                                          h₀=h_parms[3], 
#                                                          h₋₂=h_parms[1], 
#                                                          q_mult=1,
#                                                          cov_mult=1,
#                                                          R_mult=1,
#                                                          q_a=0.1,
#                                                          dynamickf=true,
#                                                          state_num=3,
#                                                          dll_b=0.1,
#                                                          acquisition_T=1e-3,
#                                                          fine_acq_T=20e-3,
#                                                          tracking_T=4e-3,
#                                                          M=20,
#                                                          σᵩ²=missing)

# figure()
# if size(trackresults.K)[1] == 3
#     Bs = get_B.(init_H.(trackresults.K[1,:], trackresults.K[2,:], 
#                         trackresults.K[3,:]))
# else
#     Bs = get_B.(init_H.(trackresults.K[1,:], trackresults.K[2,:]))
# end
# plot(trackresults.t, Bs)
# ylim([0, 10])
# xlabel("Time (s)")
# ylabel("Equivalent Bandwidth (Hz)")

# figure()
# plot(trackresults.t, trackresults.K[1,:], label="K₁")
# plot(trackresults.t, trackresults.K[2,:], label="K₂")
# if size(trackresults.K)[1] == 3
#     plot(trackresults.t, trackresults.K[3,:], label="K₃")
# end
# ylim([0,10])
# xlabel("Time (s)")
# ylabel("Kₙ (unitless)")
# legend()



# h₀ = 1e-19
# h₋₂ = 2e-20
# T = 0.001  # s
# qₐ = 0.1
# CN0 = 36 # dB⋅Hz
# sig_freq = L1_freq
# state_num = 2
# # R = [1e-4]
# R = [1/(2*10^(CN0/10)*T)]
# A = calcA(T, state_num)
# C = calcC(T, state_num)
# Q = calcQ(T, h₀, h₋₂, qₐ, sig_freq, state_num)
# K = dkalman(A, C, Q, GNSSTools.Diagonal(R))
# if length(K) == 2
#     k1 = K[1]
#     k2 = K[2]
#     B_truth = (k1 + k2/k1)/4
#     B = get_B(init_H(k1, k2))
#     println([B, B_truth, B_truth-B])
# else
#     k1 = K[1]
#     k2 = K[2]
#     k3 = K[3]
#     B = get_B(init_H(k1, k2, k3))
#     println(B)
# end

# ######################


# """
#     sim_KF_B(state_num)


# Estimates the steady state KF gain for various C/N₀ and qₐ.
# """
# function sim_KF_B(CN0s, qas, state_num; iterations=25, T=1e-3)
#     t_length = 1
#     f_s = 5e6
#     h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
#     # h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
#     sigtype = define_l1ca_code_type(; sig_freq=1*L1_freq)
#     signal = definesignal(sigtype, f_s, t_length; 
#                         prn=1, code_start_idx=1, f_d=1000, fd_rate=0, 
#                         f_if=1.25e6, CN0=36, include_phase_noise=true, 
#                         include_thermal_noise=true, include_adc=true, 
#                         include_carrier=true, receiver_h_parms=h_parms)
#     Bs = zeros(length(qas), length(CN0s))
#     p = Progress(iterations*length(CN0s)*length(qas), 1, "Processing...")
#     for j in 1:length(CN0s)
#         CN0 = CN0s[j]
#         for i in 1:length(qas)
#             for iteration in 1:iterations
#                 ϕ₀ = rand(0:0.0001:2π)
#                 f_d = rand(-5e3:0.1:5e5)
#                 # code_start_idx = rand(1:1:floor(Int, T*f_s))
#                 code_start_idx = 1
#                 definesignal!(signal; CN0=CN0, new_phase_noise=true, 
#                              new_thermal_noise=true, phi=ϕ₀, f_d=f_d,
#                              code_start_idx=code_start_idx)
#                 generatesignal!(signal)
#                 # σᵩ² = 1/(2*10^(CN0s[j]/10)*T)  # rad²
#                 σᵩ² = GNSSTools.phase_noise_variance(CN0, T)
#                 (acqresults, trackresults, corr_result, 
#                  SNR_est) = process(signal, 
#                                     sigtype, 
#                                     1, "I"; fd_range=5000, 
#                                     fine_acq_method=:fft, 
#                                     return_corrresult=true, 
#                                     h₀=h_parms[3], 
#                                     h₋₂=h_parms[1], 
#                                     q_mult=1,
#                                     σω=1000,
#                                     cov_mult=1,
#                                     R_mult=1,
#                                     q_a=qas[i],
#                                     dynamickf=true,
#                                     state_num=state_num,
#                                     dll_b=0.1,
#                                     acquisition_T=1e-3,
#                                     fine_acq_T=20e-3,
#                                     tracking_T=T,
#                                     M=20,
#                                     show_plot=false,
#                                     σᵩ²=σᵩ²)
#                 if size(trackresults.K)[1] == 3
#                     Bs[i,j] += get_B(init_H(trackresults.K[1,end], 
#                                             trackresults.K[2,end],
#                                             trackresults.K[3,end]))
#                 else
#                     Bs[i,j] += get_B(init_H(trackresults.K[1,end], 
#                                             trackresults.K[2,end]))
#                 end
#                 # Bs[i,j] += get_B(init_H(trackresults.K[:,end]...))
#                 next!(p)
#             end
#             Bs[i,j] = Bs[i,j] / iterations
#         end
#     end
#     return (Bs, CN0s, qas)
# end

# CN0s = Array(range(25, 50, step=1))
# # qas = [0.0, 0.1, 1, 10, 100]
# qas = [0, 1.0]
# Bs_2state, CN0s, qas_ = sim_KF_B(CN0s, 0.0, 2; iterations=10, T=1e-3)
# Bs_3state, CN0s, qas = sim_KF_B(CN0s, qas, 3; iterations=10, T=1e-3)


# save("ch3_kf_pll_bandwidth_3-state.jld", "Bs_2state", Bs_2state)
# save("ch3_kf_pll_bandwidth_3-state.jld", "Bs_3state", Bs_3state)
# # results = load("ch4_acquisition_performance.jld", "results")

# # figure(figsize=(7.5,2.5))
# # ax1 = subplot(1,2,1)
# # for i in 1:length(qas)
# #     ax1.plot(CN0s, Bs_2state[i,:], label="qₐ = $(qas[i])")
# # end
# # xlabel("C/N₀ (dB⋅Hz)")
# # ylabel("Bandwidth (Hz)")
# # title("(a)")
# # ax1.legend()
# figure(figsize=(7.5,3))
# ax2 = subplot(1,1,1)
# for i in reverse(1:length(qas))
#     ax2.plot(CN0s, Bs_3state[i,:], label="3-state, qₐ = $(qas[i])")
# end
# ax2.plot(CN0s, Bs_2state[1,:], "k-", label="2-state")
# xlabel("C/N₀ (dB⋅Hz)")
# ylabel("Bandwidth (Hz)")
# ax2.legend()
# subplots_adjust(wspace=0.4, bottom=0.15, left=0.08, right=0.92, top=0.95)
# savefig(string(directory, "figures/ch3_kf_pll_bandwidth.svg"), dpi=300)



################################################################################

function estimate_BKσ(T, CN0, state_num, qₐ, h_parms, sig_freq)
    h₀ = h_parms[3]
    h₋₁ = h_parms[2]
    h₋₂ = h_parms[1]
    A = calcA(T, state_num)
    C = calcC(T, state_num)
    Q = calcQ(T, h₀, h₋₂, qₐ, sig_freq, state_num)
    σᵩ²_unfiltered_expected = phase_noise_variance(CN0, T)
    R = [σᵩ²_unfiltered_expected]
    K = dkalman(A, C, Q, GNSSTools.Diagonal(R))
    # B = get_B(init_H(K...))
    B = calc_B(K...)
    σᵩ²_filtered_expected = phase_noise_variance(B, CN0, T, h₀, h₋₁, h₋₂)
    return (K, B, σᵩ²_unfiltered_expected, σᵩ²_filtered_expected)
end


function get_kf_steady_state_info(h_parms, Ts, CN0s, qₐ, state_num)
    Ks = zeros(length(h_parms), length(Ts), length(CN0s), state_num)
    Bs = zeros(length(h_parms), length(Ts), length(CN0s))
    σᵩ²s_unfiltered = zeros(length(h_parms), length(Ts), length(CN0s))
    σᵩ²s_filtered = zeros(length(h_parms), length(Ts), length(CN0s))
    p = Progress(length(h_parms)*length(CN0s)*length(Ts), 1, "Processing...")
    for i in 1:length(h_parms)
        h_parm = h_parms[i]
        for j in 1:length(Ts)
            T = Ts[j]
            for k in 1:length(CN0s)
                CN0 = CN0s[k]
                K, B, σᵩ²_unfiltered, σᵩ²_filtered = estimate_BKσ(T, CN0, 
                                                                 state_num, qₐ, 
                                                                 h_parm, 
                                                                 sig_freq)
                Ks[i,j,k,:] = K
                Bs[i,j,k] = B
                σᵩ²s_unfiltered[i,j,k] = σᵩ²_unfiltered
                σᵩ²s_filtered[i,j,k] = σᵩ²_filtered
                next!(p)
            end
        end
    end
    return (Ks, Bs, σᵩ²s_unfiltered, σᵩ²s_filtered)
end


sig_freq = L1_freq
Ts = [1e-3, 2e-3, 5e-3, 10e-3, 20e-3]  # seconds
CN0s = Array(range(20, 50, step=1))

tcxo_h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
ocxo_h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
h_parms = [tcxo_h_parms, ocxo_h_parms]


fig = figure(figsize=(7.5, 3.5))
cmap = get_cmap("viridis")

# 2-state KF
qₐ = 0
state_num = 2
Ks, Bs, σᵩ²s_unfiltered, σᵩ²s_filtered = get_kf_steady_state_info(h_parms, Ts, 
                                                                 CN0s, qₐ, 
                                                                 state_num);
# ax1 = fig.add_subplot(2,2,1)
ax1 = fig.add_subplot(1,1,1)
for i in 1:length(Ts)
    color = cmap(float(i)/length(Ts))
    ax1.plot(CN0s, Bs[1,i,:], "-", color=color,
             label="T = $(floor(Int, Ts[i]*1000))ms; TCXO")
    # ax1.plot(CN0s, Bs[2,i,:], ":", color=color,
            #  label="T = $(floor(Int, Ts[i]*1000))ms; OCXO")
end
for i in 1:length(Ts)
    color = cmap(float(i)/length(Ts))
    ax1.plot(CN0s, Bs[2,i,:], ":", color=color,
             label="T = $(floor(Int, Ts[i]*1000))ms; OCXO")
end

xlabel("C/N₀ (dB⋅Hz)")
ylabel("B (Hz)")
# title("(a)")
legend(ncol=2)

# # 3-state KF
# qₐ = 0.1
# state_num = 3
# Ks, Bs, σᵩ²s_unfiltered, σᵩ²s_filtered = get_kf_steady_state_info(h_parms, Ts, 
#                                                                  CN0s, qₐ, 
#                                                                  state_num);
# ax2 = fig.add_subplot(2,2,2)
# for i in 1:length(Ts)
#     color = cmap(float(i)/length(Ts))
#     ax2.plot(CN0s, Bs[1,i,:], "-", color=color)
#     ax2.plot(CN0s, Bs[2,i,:], ":", color=color,
#              label="T = $(floor(Int, Ts[i]*1000))ms; OCXO")
# end
# xlabel("C/N₀ (dB⋅Hz)")
# ylabel("B (Hz)")
# title("(b)")
# legend()

# # 3-state KF
# qₐ = 1
# state_num = 3
# Ks, Bs, σᵩ²s_unfiltered, σᵩ²s_filtered = get_kf_steady_state_info(h_parms, Ts, 
#                                                                  CN0s, qₐ, 
#                                                                  state_num);
# ax3 = fig.add_subplot(2,2,3)
# for i in 1:length(Ts)
#     color = cmap(float(i)/length(Ts))
#     ax3.plot(CN0s, Bs[1,i,:], "-", color=color)
#     ax3.plot(CN0s, Bs[2,i,:], ":", color=color)
# end
# xlabel("C/N₀ (dB⋅Hz)")
# ylabel("B (Hz)")
# title("(c)")


# # 3-state KF
# qₐ = 100
# state_num = 3
# Ks, Bs, σᵩ²s_unfiltered, σᵩ²s_filtered = get_kf_steady_state_info(h_parms, Ts, 
#                                                                  CN0s, qₐ, 
#                                                                  state_num);
# ax4 = fig.add_subplot(2,2,4)
# for i in 1:length(Ts)
#     color = cmap(float(i)/length(Ts))
#     ax4.plot(CN0s, Bs[1,i,:], "-", color=color)
#     ax4.plot(CN0s, Bs[2,i,:], ":", color=color)
# end
# xlabel("C/N₀ (dB⋅Hz)")
# ylabel("B (Hz)")
# title("(d)")

subplots_adjust(wspace=0.25, hspace=0.45, bottom=0.15, 
                left=0.08, right=0.92, top=0.94)

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
savefig(string(directory, "figures/ch3_kf_pll_bandwidth_numerical_2-state.svg"), dpi=300)
