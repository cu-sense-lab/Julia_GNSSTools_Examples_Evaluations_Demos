using GNSSTools
using Statistics
using ProgressMeter
using JLD
using PyPlot
pygui(true)

function get_theo_Pd(f_s, T, CN0, B, p_fa, Tsys=535)
    N = f_s * T
    v_i = N * sqrt(2*GNSSTools.k*535)*10^(CN0/20)
    σ_n_theoretical = sqrt(GNSSTools.k*B*Tsys) 
    SNR_theo = CN0 - 10*log10(1/T)
    v_t_theoretical = p_fa2v_t(p_fa, sqrt(N)*σ_n_theoretical)
    P_d_theo = v_t2p_d(v_t_theoretical, 
                       sqrt(N)*σ_n_theoretical, 
                       v_i=v_i)
    # P_d_theo = v_t2p_d(v_t_theoretical, sqrt(N)*σ_n_theoretical, 
                    #    SNR=SNR_theo)
    return P_d_theo
end

"""
sim_pd(CN0s, state_num; iterations=100, T=2e-3)


Simulates signals and processes them to calculate the propability of detection,
P_d.
"""
function sim_pd(CN0s, f_s; iterations=100, T=1e-3, M=1, 
                return_theoretical=true, f_d=0,
                sigtype=define_l1ca_code_type(),
                doppler_bin_num=5)
    t_length = T*M
    prn = 1
    fd_range = 5e3
    # f_if = 1.25e6
    f_if = 0
    h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
    # h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
    # sigtype = define_l1ca_code_type()
    # sigtype = define_l5_code_type()
    if sigtype.include_I && sigtype.include_Q
        B = max(sigtype.B_I, sigtype.B_Q)
    elseif sigtype.include_I && !sigtype.include_Q
        B = sigtype.B_I
    elseif !sigtype.include_I && sigtype.include_Q
        B = sigtype.B_Q
    else
        error("No CodeType defined in SigType struct.")
    end
    d = ceil(Int, f_s/l1ca_chipping_rate/2)
    replica = definereplica(sigtype, f_s, T)
    signal = definesignal(sigtype, f_s, t_length; 
                          prn=prn, fd_rate=0, f_if=f_if,
                          include_phase_noise=false, include_thermal_noise=true, 
                          include_adc=true, include_carrier=true, 
                          receiver_h_parms=h_parms)
    if return_theoretical
        P_ds_theoretical = Array{Float64}(undef, length(CN0s))
        SNR_theoretical = Array{Float64}(undef, length(CN0s))
    end
    P_ds_experimental = Array{Float64}(undef, length(CN0s))
    P_ds_exp_from_dopler_code = Array{Float64}(undef, length(CN0s))
    SNR_experimental = Array{Float64}(undef, length(CN0s))
    p = Progress(iterations*length(CN0s), 1, "Processing...")
    for j in 1:length(CN0s)
        σ²_ns = zeros(iterations)
        powers = zeros(iterations)
        above_threshold_num_1 = 0
        above_threshold_num_2 = 0
        SNR_avg = 0.
        for iteration in 1:iterations
            # ϕ₀ = rand(0:0.0001:2π)
            ϕ₀ = 0
            # f_d = rand(-fd_range:0.1:fd_range)
            # code_start_idx = rand(1:1:replica.sample_num)
            code_start_idx = 1
            # ϕ₀ = 0
            # code_start_idx = 1000
            definesignal!(signal; CN0=CN0s[j], new_phase_noise=false, 
                          new_thermal_noise=true, phi=ϕ₀, f_d=f_d,
                          code_start_idx=code_start_idx,
                          nADC=16)
            generatesignal!(signal)
            fd_range = (doppler_bin_num - 1)/(2*T)
            fd_course, n0_est, SNR_est, P_d, above_threshold,
            corr_result = courseacquisition(signal, 
                                            replica,
                                            prn; 
                                            fd_center=0, 
                                            fd_range=fd_range,
                                            fd_rate=0,
                                            return_corrresult=true,
                                            M=M)
            rows, cols = size(corr_result)
            idx = argmax(corr_result)
            row_idx = idx[1]
            col_idx = idx[2]
            noise_row_idx = (row_idx + 3) % rows + 1
            # Peak power should be peak + the two adjacent bins
            # lower_row_idx = row_idx - 1
            # if lower_row_idx < 1
            #     lower_row_idx = rows + lower_row_idx
            # end
            # upper_row_idx = row_idx + 1
            # if upper_row_idx > rows
            #     upper_row_idx = upper_row_idx % rows
            # end
            # lower_col_idx = col_idx - 1
            # if lower_col_idx < 1
            #     lower_col_idx = replica.sample_num + lower_col_idx
            # end
            # upper_col_idx = col_idx + 1
            # if upper_col_idx > replica.sample_num 
            #     upper_col_idx = upper_col_idx % replica.sample_num
            # end
            # power = sum(sqrt.([corr_result[row_idx,lower_col_idx],
            #                   corr_result[row_idx,col_idx], 
            #                   corr_result[row_idx,upper_col_idx]]))
            # power = sum(sqrt.([corr_result[lower_row_idx,col_idx],
            #                   corr_result[row_idx,col_idx], 
            #                   corr_result[upper_row_idx,col_idx]]))
            power = sqrt(corr_result[argmax(corr_result)])
            σ²_n = sum(corr_result[noise_row_idx,:]) / (cols - 1)
            σ_n = sqrt(σ²_n)
            # σ²_n = (sum(corr_result) - maximum(corr_result))/(rows*cols - 1)
            # v_t = p_fa2v_t(p_fa, sqrt(σ²_n))
            v_t = p_fa2v_t(p_fa, σ_n)
            # SNR_avg += power/σ_n
            powers[iteration] = power
            σ²_ns[iteration] = σ_n
            if power >= v_t
                above_threshold_num_2 += 1
            end
            correct_doppler = fd_course == (round(f_d*T, digits=0)/T)
            # correct_n0 = (n0_est >= code_start_idx-d) && (n0_est <= code_start_idx+d)
            correct_n0 = n0_est == code_start_idx
            # correct_n0 = n0_est == code_start_idx
            if correct_doppler && correct_n0
                above_threshold_num_1 += 1
            end
            # println([n0_est, code_start_idx, f_d, fd_course])
            next!(p)
        end
        # σ_n_exp = mean(σ²_ns)
        # # println([minimum(σ²_ns), maximum(σ²_ns), maximum(σ²_ns)-minimum(σ²_ns)])
        # v_t_exp = p_fa2v_t(p_fa, σ_n_exp)
        # for iteration in 1:iterations
        #     if powers[iteration] >= v_t_exp
        #         above_threshold_num_2 += 1
        #     end
        # end

        SNR_experimental[j] = 10*log10(SNR_avg / iterations)
        P_ds_experimental[j] = above_threshold_num_2 / iterations
        P_ds_exp_from_dopler_code[j] = above_threshold_num_1 nbins = 1000
        fig = figure(figsize=(7.5, 4.5))
        ax1 = fig.add_subplot(2,2,1, aspect="auto")
        hist2D(gps_hist.dopplers./1000, gps_hist.doppler_rates, bins=nbins)
        xlabel("Doppler (kHz)")
        ylabel("Doppler Rate (Hz/s)")
        title("(a)")
        
        ax2 = fig.add_subplot(2,2,2, aspect="auto")
        hist2D(iridium_hist.dopplers./1000, iridium_hist.doppler_rates, bins=nbins)
        xlabel("Doppler (kHz)")
        ylabel("Doppler Rate (Hz/s)")
        title("(b)")
        
        ax3 = fig.add_subplot(2,2,3, aspect="auto")
        hist2D(starlink_hist.dopplers./1000, starlink_hist.doppler_rates, bins=nbins)
        xlabel("Doppler (kHz)")
        ylabel("Doppler Rate (Hz/s)")
        title("(c)")
        
        ax4 = fig.add_subplot(2,2,4, aspect="auto")
        hist2D(oneweb_hist.dopplers./1000, oneweb_hist.doppler_rates, bins=nbins)
        xlabel("Doppler (kHz)")
        ylabel("Doppler Rate (Hz/s)")
        title("(d)")/ iterations
        if return_theoretical
            P_d_theo = get_theo_Pd(f_s, T, CN0s[j], B, p_fa, signal.Tsys)
            # σ_n_theoretical = sqrt(GNSSTools.k*B*signal.Tsys)
            # # v_i = f_s * T * (sqrt(2*GNSSTools.k*535)*10^(CN0s[j]/20))^2
            # N = f_s*T
            # v_i = N * sqrt(GNSSTools.k*535)*10^(CN0s[j]/20)
            # SNR_theo = CN0s[j] - 10*log10(1/T)
            # SNR_theoretical[j] = SNR_theo
            # v_t_theoretical = p_fa2v_t(p_fa, sqrt(N)*σ_n_theoretical)
            # P_ds_theoretical[j] = v_t2p_d(v_t_theoretical, 
            #                               sqrt(N)*σ_n_theoretical, 
            #                               v_i=v_i)
            # # P_ds_theoretical[j] = v_t2p_d(v_t_theoretical, σ_n_theoretical, 
            #                             #   SNR=SNR_theo)
            P_ds_theoretical[j] = P_d_theo
        end
    end
    if return_theoretical
        return (P_ds_theoretical, P_ds_experimental, P_ds_exp_from_dopler_code)
    else
        return (P_ds_experimental, P_ds_exp_from_dopler_code)
    end
end

CN0s = Array(range(25, 50, step=0.5))
iterations = 100
p_fa = 1e-7
f_d = 0

# f_s = 2.5e6
# sigtype = define_l1ca_code_type()
# Ts = [1e-3 2e-3 5e-3 10e-3]
# name = "l1ca"

f_s = 25e6
sigtype = define_l5_code_type(;channel="Q")
Ts = [1e-3, 20e-3]  # seconds
name = "l5"

M = 1 
P_ds_theoretical = []
P_ds_experimental = []
P_ds_exp_from_dopler_code = []

println("Working on $(sigtype.name) signal...")

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name = "ch4_compare_theoretical_empirical_Pd_$(name)_fd=$(f_d)Hz"
file = string(directory, file_name, ".jld")

for i in 1:length(Ts)
        P_ds_theo, P_ds_exp, P_ds_exp2 = sim_pd(CN0s, f_s; iterations=iterations, 
                                                T=Ts[i], M=M, 
                                                return_theoretical=true,
                                                sigtype=sigtype,
                                                f_d=f_d,
                                                doppler_bin_num=5)
        push!(P_ds_theoretical, P_ds_theo)
        push!(P_ds_experimental, P_ds_exp)
        push!(P_ds_exp_from_dopler_code, P_ds_exp2)
        save(file, 
             "P_ds_theoretical", P_ds_theoretical, 
             "P_ds_experimental", P_ds_experimental,
             "P_ds_exp_from_dopler_code", P_ds_exp_from_dopler_code,
             "CN0s", CN0s,
             "M", M,
             "p_fa", p_fa,
             "iterations", iterations,
             "Ts", Ts,
             "f_s", f_s)
end

P_ds_theoretical, P_ds_experimental, CN0s, f_s = load(file, 
                                                 "P_ds_theoretical", 
                                                 "P_ds_experimental",
                                                 "CN0s",
                                                 "f_s")

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_names = [string(directory, "ch4_compare_theoretical_empirical_Pd_l1ca_fd=0Hz.jld"),
              string(directory, "ch4_compare_theoretical_empirical_Pd_l5_fd=0Hz.jld")]
cmap = get_cmap("viridis")
fig = figure(figsize=(7.5, 3.5))
for f in 1:length(file_names)
    file_name = file_names[f]
    P_ds_theoretical, P_ds_experimental, CN0s, f_s, Ts = load(file_name, 
                                                          "P_ds_theoretical", 
                                                          "P_ds_experimental",
                                                          "CN0s",
                                                          "f_s",
                                                          "Ts")
    Pds_to_plot = length(P_ds_theoretical)
    Pd_to_plot = Array(1:Pds_to_plot)
    ax = fig.add_subplot(1, 2, f)
    for i in 1:Pds_to_plot
        if i in Pd_to_plot
            color = cmap(float(i)/Pds_to_plot)
            label = "T = $(floor(Int, Ts[i]*1e3))ms"
            # label = "T = $(floor(Int, Ts[i]*1e3))ms, f_s = $(round(f_s*1e-6, digits=1))MHz"
            ax.plot(CN0s, P_ds_theoretical[i], "-", color=color, 
                    label=label)
            ax.plot(CN0s, P_ds_experimental[i], ":", color=color)
            ax.plot(CN0s, P_ds_experimental[i], ".", color=color)
        end
    end
    legend()
    xlabel("C/N₀ (dB⋅Hz)")
    ylabel(L"P_d")
    letter = ('a' + (f-1))
    title("($(letter))")
end
subplots_adjust(wspace=0.3, bottom=0.15, left=0.08, right=0.92, top=0.93)
savefig(string(directory, "figures/ch4_compare_theoretical_empirical_Pd.svg"), dpi=300)

# """

# The error between the truth and observation is near zero when f_d is an integer 
# multiple of 1/T. If f_d is not a multiple of 1/T and is 1/2T from either 
# close multiples of 1/T, the error will be the maximum (~3.5dB).

# Also, if f_s is a multiple of the code chipping rate, then it will also
# result in the worst acquisition performance and greater error between the
# the truth and observation.

# """