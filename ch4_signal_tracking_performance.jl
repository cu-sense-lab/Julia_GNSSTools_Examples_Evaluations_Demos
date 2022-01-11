using GNSSTools
using Statistics
using ProgressMeter
using JLD
using LinearAlgebra
using Optim
using Distributions
using QuadGK
using PyPlot
pygui(true)
include("kalman_filter_gain_error_functions.jl")


function eval_tracking(CN0s, f_s; 
                       iterations=iterations, 
                       T=1e-3, M=5, p_fa=1e-7,
                       state_num=3,
                       sigtype=define_l1ca_code_type(),
                       channel="I",
                       h_parms=[1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2,
                       qₐ=1e-9,
                       σω=1,
                       acquisition_T=20e-3,
                       fine_acq_T=100e-3,
                       t_length=1,
                       dll_b=1,
                       fd_rate=0,
                       f_if=0,
                       P0=missing)
    prn = 1
    fd_range = 2/acquisition_T
    sig_freq = sigtype.sig_freq
    h₀ = h_parms[3]
    h₋₁ = h_parms[2]
    h₋₂ = h_parms[1]
    replica = definereplica(sigtype, f_s, T)
    signal = definesignal(sigtype, f_s, t_length; 
                          prn=prn, fd_rate=fd_rate, f_if=f_if,
                          include_phase_noise=true, 
                          include_thermal_noise=true, 
                          include_adc=true, include_carrier=true, 
                          receiver_h_parms=h_parms)
    σᵩ_theo = Array{Float64}(undef, length(CN0s))
    σᵩ_exp = Array{Float64}(undef, length(CN0s))
    σᵩ_exp_std = Array{Float64}(undef, length(CN0s))
    σ_unfilt = Array{Float64}(undef, length(CN0s))
    σ_unfilt_std = Array{Float64}(undef, length(CN0s))
    pᵩs = Array{Float64}(undef, length(CN0s))
    σᵩ²_unfilt_theo = Array{Float64}(undef, length(CN0s))
    pᵩs_steady_state = Array{Float64}(undef, length(CN0s))
    Bs = Array{Float64}(undef, length(CN0s))
    Ps = Array{Float64}(undef, state_num, length(CN0s))
    Ks = Array{Float64}(undef, state_num, length(CN0s))
    P_fulls = Array{Float64}(undef, length(CN0s), state_num, state_num)
    N_over_2 = floor(Int, t_length/T/2)
    p = Progress(iterations*length(CN0s), 1, "Processing...")
    for i in 1:length(CN0s)
        CN0 = CN0s[i]
        CN0_linear = 10^(CN0/10)
        σᵩ² = GNSSTools.phase_noise_variance(CN0, T)
        R = [σᵩ²]
        if ismissing(P0)
            if state_num == 3
                P0 = diagm([σᵩ², (2π/(2*T))^2, (2π*σω)^2])
            elseif state_num == 2
                P0 = diagm([σᵩ², (2π/(2*T))^2])
            else
                @error("Invalid value for state_num.")
            end
        end
        sigma_phis = Array{Float64}(undef, iterations)
        pᵩ = Array{Float64}(undef, iterations)
        P = Array{Float64}(undef, state_num, iterations)
        P_full = Array{Float64}(undef, iterations, state_num, state_num)
        K = Array{Float64}(undef, state_num, iterations)
        sigma_unfiltered = Array{Float64}(undef, iterations)
        for iteration in 1:iterations
            ϕ₀ = rand(0:0.0001:2π)
            f_d = rand(-fd_range:0.1:fd_range)
            if state_num == 3
                x = [ϕ₀, 2π*f_d, 2π*fd_rate] #.+ rand(MvNormal(P0))
                ϕ0_est = x[1]
                fd_est = x[2]/2π
                fd_rate_est = x[3]/2π
            else
                x = [ϕ₀, 2π*f_d] #.+ rand(MvNormal(P0))
                ϕ0_est = x[1]
                fd_est = x[2]/2π
            end
            # code_start_idx = rationalize(rand(1:0.25:replica.sample_num))
            code_start_idx = rand(1:1:replica.sample_num)
            # code_start_idx = 1
            definesignal!(signal; CN0=CN0, new_phase_noise=true, 
                          new_thermal_noise=true, phi=ϕ₀, f_d=f_d,
                          code_start_idx=code_start_idx,
                          nADC=16, fd_rate=fd_rate)
            generatesignal!(signal)
            replica = definereplica(sigtype, f_s, T)
            trackresults = trackprn(signal, replica, prn, ϕ0_est,
                                    fd_est, 
                                    code_start_idx, 
                                    P0, R; 
                                    DLL_B=dll_b,
                                    state_num=state_num, 
                                    dynamickf=false,
                                    cov_mult=1, 
                                    qₐ=qₐ, 
                                    q_mult=1,
                                    h₀=h₀, h₋₂=h₋₂, 
                                    R_mult=1, 
                                    fd_rate=0)
            # acqresults, trackresults, corr_result, SNR_est = process(signal, 
            #                             sigtype, 
            #                             prn, channel;
            #                             fine_acq_method=:fft, 
            #                             return_corrresult=true, 
            #                             fd_center=f_d,
            #                             fd_range=fd_range, 
            #                             h₀=h₀, 
            #                             h₋₂=h₋₂, 
            #                             q_mult=1,
            #                             cov_mult=1,
            #                             R_mult=1,
            #                             σᵩ²=σᵩ²,
            #                             σω=σω,
            #                             q_a=qₐ,
            #                             dynamickf=true,
            #                             state_num=state_num,
            #                             dll_b=dll_b,
            #                             acquisition_T=acquisition_T,
            #                             fine_acq_T=fine_acq_T,
            #                             tracking_T=T,
            #                             M=M,
            #                             show_plot=false)
            K[:,iteration] = mean(trackresults.K[:,N_over_2:end], dims=2)
            sigma_phis[iteration] = var(trackresults.dphi_filt[N_over_2:end])
            sigma_unfiltered[iteration] = var(trackresults.dphi_meas[N_over_2:end])
            P[:,iteration] = mean(trackresults.P[:,N_over_2:end], dims=2)
            P_full[iteration,:,:] = mean(trackresults.P_full[N_over_2:end,:,:], dims=1)
            C = calcC(T, state_num)
            pᵩ[iteration] = (C*P_full[iteration,:,:]*C')[1]
            next!(p)
        end
        pᵩs[i] = mean(pᵩ)
        σᵩ_exp[i] = sqrt(mean(sigma_phis))
        σᵩ_exp_std[i] = std(sqrt.(sigma_phis))
        σ_unfilt[i] = sqrt(mean(sigma_unfiltered))
        σ_unfilt_std[i] = std(sqrt.(sigma_unfiltered))
        Ps[:,i] = mean(P, dims=2)
        P_fulls[i,:,:] = mean(P_full, dims=1)
        K_mean = mean(K, dims=2)
        Ks[:,i] = K_mean
        B = get_B(init_H(K_mean...))
        A = calcA(T, state_num)
        C = calcC(T, state_num)
        Q = calcQ(T, h₀, h₋₂, qₐ, sig_freq, state_num)
        R = [σᵩ²]
        σᵩ²_unfilt_theo[i] = σᵩ²
        S = GNSSTools.dare(A', C', Q, GNSSTools.Diagonal(R)')
        # P_steady_state = A*S*A'
        P_steady_state = S
        pᵩ_steady_state = C*P_steady_state*C'
        pᵩs_steady_state[i] = pᵩ_steady_state[1]
        Bs[i] = B
        σᵩ_theo[i] = sqrt(phase_noise_variance(B, CN0, T, h₀, h₋₁, h₋₂))
    end
    return (σᵩ_theo, σᵩ_exp, Ps, Bs, σ_unfilt, Ks, pᵩs, P_fulls, pᵩs_steady_state,
            σᵩ²_unfilt_theo, σᵩ_exp_std, σ_unfilt_std)
end

h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
# h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
h₀ = h_parms[3]
h₋₁ = h_parms[2]
h₋₂ = h_parms[1]

channel = "I"
sigtype = define_l1ca_code_type()
# Ts = [1e-3 2e-3 5e-3 10e-3]
Ts = [1e-3 10e-3]
f_s = 2.5e6
name = "l1ca"

# channel = "Q"
# sigtype = define_l5_code_type(;channel=channel)
# Ts = [1e-3 20e-3]
# f_s = 25e6
# name = "l5"

CN0s = Array(range(20, 50, step=1))
iterations = 25
p_fa = 1e-7
qₐ = 0.1
state_num = 2
σω = 1
acquisition_T=40e-3
fine_acq_T=100e-3
M=2
t_length=1
dll_b=1
fd_rate = 0

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name = "ch4_signal_tracking_performance_$(name)"
file = string(directory, file_name, ".jld")

σ_theoretical = []
σ_experimental = []
σ_unfiltered = []
Ks = []
Ps = []
pᵩs = []
Bs = []
P_fulls = []
pᵩs_steady_state = []
σᵩ²_unfilt_theo = []
σᵩ_exp_std = []
σ_unfilt_std = []

for i in 1:length(Ts)
    σ_theo, σ_exp, P, B, σ_unfilt, K, pᵩ, P_full, 
    pᵩ_steady_state, σᵩ²_unfilt_theo_, σᵩ_exp_std_, σ_unfilt_std_  = eval_tracking(CN0s, f_s; 
                                                     iterations=iterations, 
                                                     T=Ts[i], p_fa=p_fa,
                                                     state_num=state_num,
                                                     sigtype=sigtype,
                                                     channel=channel,
                                                     qₐ=qₐ,
                                                     σω=σω,
                                                     acquisition_T=acquisition_T,
                                                     fine_acq_T=fine_acq_T,
                                                     M=M,
                                                     t_length=t_length,
                                                     dll_b=dll_b,
                                                     fd_rate=fd_rate)
    push!(σ_theoretical, σ_theo)
    push!(σ_experimental, σ_exp)
    push!(Ps, P)
    push!(Bs, B)
    push!(Ks, K)
    push!(σ_unfiltered, σ_unfilt)
    push!(pᵩs, pᵩ)
    push!(P_fulls, P_full)
    push!(pᵩs_steady_state, pᵩ_steady_state)
    push!(σᵩ²_unfilt_theo, σᵩ²_unfilt_theo_)
    push!(σᵩ_exp_std, σᵩ_exp_std_)
    push!(σ_unfilt_std, σ_unfilt_std_)
    save(file, 
         "σ_theoretical", σ_theoretical,
         "σ_experimental", σ_experimental,
         "Ps", Ps,
         "Bs", Bs,
         "Ts", Ts,
         "CN0s", CN0s,
         "iterations", iterations,
         "p_fa", p_fa,
         "f_s", f_s,
         "M", M,
         "h_parms", h_parms,
         "sig_freq", sigtype.sig_freq,
         "q_a", qₐ,
         "state_num", state_num,
         "σω", σω,
         "acquisition_T", acquisition_T,
         "fine_acq_T", fine_acq_T,
         "t_length", t_length,
         "dll_b", dll_b,
         "Ks", Ks,
         "σ_unfiltered", σ_unfiltered,
         "pᵩs", pᵩs,
         "P_fulls", P_fulls,
         "pᵩs_steady_state", pᵩs_steady_state,
         "fd_rate", fd_rate,
         "σᵩ²_unfilt_theo", σᵩ²_unfilt_theo,
         "σᵩ_exp_std", σᵩ_exp_std,
         "σ_unfilt_std", σ_unfilt_std)
end


directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_names = [string(directory, "ch4_signal_tracking_performance_l1ca.jld"),
              string(directory, "ch4_signal_tracking_performance_l5.jld")]
cmap = get_cmap("viridis")

fig = figure(figsize=(7.5, 3))
for f in 1:1
# for f in 1:length(file_names)
    file_name = file_names[f]
    ax = fig.add_subplot(1, 2, f)
    CN0s,σ_theoretical,σ_experimental, Ts, σ_unfiltered, σᵩ_exp_std, σ_unfilt_std = load(file_name, 
                                             "CN0s", 
                                             "σ_theoretical",
                                             "σ_experimental",
                                             "Ts",
                                             "σ_unfiltered",
                                             "σᵩ_exp_std",
                                             "σ_unfilt_std")
    σs_to_plot = length(σ_theoretical)
    σ_to_plot = Array(1:σs_to_plot)
    for i in 1:σs_to_plot
        if i in σ_to_plot
            color = cmap(float(i)/σs_to_plot)
            # label = "T = $(floor(Int, Ts[i]*1e3))ms, f_s = $(round(f_s*1e-6, digits=1))MHz"
            label = "T = $(floor(Int, Ts[i]*1e3))ms"
            # ax.plot(CN0s, sqrt.(pᵩs_steady_state[i]).*(180/π), "-", color=color, 
                    # label=label)
            # ax.plot(CN0s, σ_theoretical[i].*(180/π), ":", color=color)
            # ax.plot(CN0s, σ_experimental[i].*(180/π), ":", color=color)
            # ax.plot(CN0s, σ_experimental[i].*(180/π), "+", color=color)
            # ax.plot(CN0s, sqrt.(σᵩ²_unfilt_theo[i]).*(180/π), "-", color=color)
            # ax.plot(CN0s, sqrt.(Ps[i][1,:]).*(180/π), ":", color=color)
            # ax.plot(CN0s, sqrt.(Ps[i][1,:]).*(180/π), ".", color=color)
            # ax.plot(CN0s, sqrt.(pᵩs[i]).*(180/π), ".", color=color)
            ax.plot(CN0s, sqrt.(GNSSTools.phase_noise_variance.(CN0s, Ts[i])).*(180/π), 
                    "-", color=color)
            ax.errorbar(CN0s, σ_unfiltered[i].*(180/π), yerr=σ_unfilt_std[i],
                        marker=".", color=color)
            ax.plot([], [], ".-", color=color, 
                    label="T=$(floor(Int, Ts[i]*1000))ms")
        end
    end
    ax.set_yscale("log")
    # hlines(y=10log10(15), xmin=CN0s[1], xmax=CN0s[end], color="grey", 
        #    linestyle=":", label=string(L"15^\circ", " Tracking Limit"))
    legend()
    xlabel("C/N₀ (dB⋅Hz)")
    ylim([0.0, 60])
    ylabel("Unfiltered σᵩ (degrees)")
    letter = ('a' + (f-1))
    title("($(letter))")
end
subplots_adjust(wspace=0.3, bottom=0.15, left=0.08, right=0.92, top=0.92)
# savefig(string(directory, "figures/ch4_signal_tracking_performance_unfiltered_sigma-phi.pdf"), dpi=300)


fig = figure(figsize=(7.5, 3))
for f in 1:1
# for f in 1:length(file_names)
    file_name = file_names[f]
    ax = fig.add_subplot(1, 2, f)
    CN0s,σ_theoretical,σ_experimental, Ts, σ_unfiltered, 
    σᵩ_exp_std, σ_unfilt_std, pᵩs_steady_state = load(file_name, 
                                             "CN0s", 
                                             "σ_theoretical",
                                             "σ_experimental",
                                             "Ts",
                                             "σ_unfiltered",
                                             "σᵩ_exp_std",
                                             "σ_unfilt_std",
                                             "pᵩs_steady_state")
    σs_to_plot = length(σ_theoretical)
    σ_to_plot = Array(1:σs_to_plot)
    for i in 1:σs_to_plot
        if i in σ_to_plot
            color = cmap(float(i)/σs_to_plot)
            # label = "T = $(floor(Int, Ts[i]*1e3))ms, f_s = $(round(f_s*1e-6, digits=1))MHz"
            label = "T = $(floor(Int, Ts[i]*1e3))ms"
            # ax.plot(CN0s, sqrt.(pᵩs_steady_state[i]).*(180/π), "-", color=color, 
                    # label=label)
            ax.plot(CN0s, sqrt.(pᵩs_steady_state[i]).*(180/π), "-", color=color)
            # ax.plot(CN0s, sqrt.(GNSSTools.phase_noise_variance.(CN0s, Ts[i])).*(180/π), ":", color=color)
            # ax.plot(CN0s, σ_experimental[i].*(180/π), ":", color=color)
            ax.errorbar(CN0s, σ_experimental[i].*(180/π), yerr=σᵩ_exp_std[i], 
                        marker=".", color=color)
            # ax.plot(CN0s, σ_unfiltered[i].*(180/π), "x", color=color)
            # ax.plot(CN0s, sqrt.(σᵩ²_unfilt_theo[i]).*(180/π), "-", color=color)
            # ax.plot(CN0s, sqrt.(Ps[i][1,:]).*(180/π), ":", color=color)
            # ax.plot(CN0s, sqrt.(Ps[i][1,:]).*(180/π), ".", color=color)
            # ax.plot(CN0s, sqrt.(pᵩs[i]).*(180/π), ".", color=color)
            ax.plot([], [], ".-", color=color, 
                    label="T=$(floor(Int, Ts[i]*1000))ms")
        end
    end
    ax.set_yscale("log")
    hlines(y=10log10(15), xmin=CN0s[1], xmax=CN0s[end], color="grey", 
           linestyle=":", label=string(L"15^\circ", " Tracking Limit"))
    legend()
    xlabel("C/N₀ (dB⋅Hz)")
    ylim([0.0, 60])
    ylabel("Filtered σᵩ (degrees)")
    letter = ('a' + (f-1))
    title("($(letter))")
end
subplots_adjust(wspace=0.3, bottom=0.15, left=0.08, right=0.92, top=0.92)
# savefig(string(directory, "figures/ch4_signal_tracking_performance_filtered_sigma-phi.pdf"), dpi=300)
