using GNSSTools
using Statistics
using Random
using StatsBase
using LinearAlgebra
using ProgressMeter
using Base.Threads
using JLD
using PyPlot
pygui(true)
include("ch4_constellation_hist_structs.jl")


function benchmark_acquisition(sigtype::GNSSTools.SignalType, Ts, CN0s, 
                               fd_rates, f_s; iterations=100,
                               include_phase_noise=true, 
                               include_thermal_noise=true,
                               include_adc=true, 
                               include_carrier=true, 
                               h_parms=[1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2,
                               code_start_idx=missing, phi=missing, f_d=missing,
                               nADC=8, fd_ranges=missing, f_if=1.25e6, prn=1, M=1,
                               Pd_threshold=0.9, p_fa=1e-7,
                               doppler_bin_num=5) 
    if sigtype.include_I && sigtype.include_Q
        B = max(sigtype.B_I, sigtype.B_Q)
    elseif sigtype.include_I && !sigtype.include_Q
        B = sigtype.B_I
    elseif !sigtype.include_I && sigtype.include_Q
        B = sigtype.B_Q
    else
        error("No CodeType defined in SigType struct.")
    end
    Pds = Array{Float64}(undef, length(Ts), length(CN0s), length(fd_rates))
    p = Progress(length(Ts)*length(CN0s)*length(fd_rates)*iterations, 1,
                 "Processing")
    @threads for i in 1:length(Ts)
        T = Ts[i]
        t_length = T*M
        replica = definereplica(sigtype, f_s, T)
        signal = definesignal(sigtype, f_s, t_length; 
                              prn=prn, f_if=f_if, nADC=nADC,
                              include_phase_noise=include_phase_noise, 
                              include_thermal_noise=include_thermal_noise, 
                              include_adc=include_adc, 
                              include_carrier=include_carrier, 
                              receiver_h_parms=h_parms)
        for j in 1:length(CN0s)
            CN0 = CN0s[j]
            for k in 1:length(fd_rates)
                fd_range = fd_ranges[k]
                fd_rate = fd_rates[k]
                Pd = 0
                for iteration in 1:iterations
                    if ismissing(phi)
                        ϕ₀_truth = rand(0:0.0001:2π)
                    else
                        ϕ₀_truth = phi
                    end
                    if ismissing(f_d)
                        f_d_truth = rand(-fd_range:0.1:fd_range)
                    else
                        f_d_truth = f_d
                    end
                    if ismissing(code_start_idx)
                        n0_truth = rand(1:1:replica.sample_num)
                    else
                        n0_truth = code_start_idx
                    end
                    if include_phase_noise
                        new_phase_noise = true
                    else
                        new_phase_noise = false
                    end
                    if include_thermal_noise
                        new_thermal_noise = true
                    else
                        new_thermal_noise = false
                    end
                    definesignal!(signal; CN0=CN0, phi=ϕ₀_truth, f_d=f_d_truth,
                                  code_start_idx=n0_truth, fd_rate=fd_rate,
                                  new_phase_noise=new_phase_noise, 
                                  new_thermal_noise=new_thermal_noise)
                    generatesignal!(signal)
                    acq_fd_range = (doppler_bin_num-1)/(2*T)
                    fd_center = round(f_d_truth*T, digits=0)/T
                    fd_course, n0_est, SNR_est, Pd_est, above_threshold,
                    corr_result = courseacquisition(signal, 
                                                    replica,
                                                    prn; 
                                                    fd_center=fd_center, 
                                                    fd_range=acq_fd_range,
                                                    fd_rate=0,
                                                    return_corrresult=true,
                                                    M=M,
                                                    p_fa=p_fa)
                    if above_threshold
                        Pd += 1
                    end
                    next!(p)
                end
                Pds[i,j,k] = Pd / iterations
            end
        end
    end
    return Pds
end


function plot_data_for_T(i, fig, Ts, CN0s, fd_rates, Pds; 
                         cmap=get_cmap("viridis"), letter=('a' + (i-1)),
                         plot_title="($(letter))")
    ax = fig.add_subplot(2, 2, i)
    for j in 1:length(Ts)
        T = Ts[j]
        color = cmap(float(j)/length(Ts))
        # 0Hz/s Doppler rate
        if i == 1
            label1 = string("T = $(floor(Int, T*1e3))ms, Low ", L"\dot{f}_d")
            ax.plot(CN0s, Pds[j,:,1], color=color, label=label1, linestyle="-")
        else
            ax.plot(CN0s, Pds[j,:,1], color=color, linestyle="-")
        end
        # Largest minimum Doppler rate
        if i == 2
            label2 = string("T = $(floor(Int, T*1e3))ms, High ", L"\dot{f}_d")
            ax.plot(CN0s, Pds[j,:,2], color=color, label=label2, linestyle=":")
        else
            ax.plot(CN0s, Pds[j,:,2], color=color, linestyle=":")
        end
    end
    xlabel("C/N₀ (dB⋅Hz)")
    ylabel(L"P_d")
    title(plot_title)
    return ax
end

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
directory = ""
doppler_file_name = string(directory, "ch1_doppler_curves.jld")
doppler_raw_gps, doppler_rate_raw_gps, sig_freq = load(doppler_file_name, 
                                                       "doppler_raw_gps",
                                                       "doppler_rate_raw_gps",
                                                       "freq")
doppler_raw_iridium, doppler_rate_raw_iridium = load(doppler_file_name, 
                                                     "doppler_raw_iridium",
                                                     "doppler_rate_raw_iridium")
doppler_raw_starlink, doppler_rate_raw_starlink = load(doppler_file_name, 
                                                       "doppler_raw_starlink",
                                                       "doppler_rate_raw_starlink")
doppler_raw_oneweb, doppler_rate_raw_oneweb = load(doppler_file_name, 
                                                   "doppler_raw_oneweb",
                                                   "doppler_rate_raw_oneweb")


gps_hist = initiate_hist(doppler_raw_gps, doppler_rate_raw_gps; nbins=1000)
iridium_hist = initiate_hist(doppler_raw_iridium, doppler_rate_raw_iridium; 
                             nbins=1000)
starlink_hist = initiate_hist(doppler_raw_starlink, doppler_rate_raw_starlink; 
                              nbins=1000)
oneweb_hist = initiate_hist(doppler_raw_oneweb, doppler_rate_raw_oneweb; 
                            nbins=1000)


iterations = 100
file_names = ["ch4_section3_acquisition_l1ca_constellation_comparison.jld",
              "ch4_section3_acquisition_l5q_constellation_comparison.jld"]
f_ss = [2.5e6, 25e6]  # Hz
f_ifs = [0e6, 0e6]  # Hz
phis = [missing, missing]
code_start_idxs = [missing, missing]
Tss = [[1, 2, 5, 10], 
       [1, 20, 40, 80, 100]] .* 1e-3  # ms
CN0s = Array(range(25, 50, step=1))
sigtypes = [define_l1ca_code_type(maximum(Tss[1])), 
            define_l5_code_type(;channel="Q")]

# items = [1, 2]
items = [2]
for i in items
    println("Will compare constellations for $(sigtypes[i].name) signal.")
end
for i in items
    f_s = f_ss[i]
    f_if = f_ifs[i]
    phi = phis[i]
    Ts = Tss[i]
    code_start_idx = code_start_idxs[i]
    sigtype = sigtypes[i]
    fd_ranges = [maximum(gps_hist.f_d),      maximum(gps_hist.f_d),
                 maximum(iridium_hist.f_d),  maximum(iridium_hist.f_d),
                 maximum(starlink_hist.f_d), maximum(starlink_hist.f_d),
                 maximum(oneweb_hist.f_d),   maximum(oneweb_hist.f_d)] .* (sigtype.sig_freq/sig_freq)
    fd_rates = [0, minimum(gps_hist.fd_rate), 
                0, minimum(iridium_hist.fd_rate),
                0, minimum(starlink_hist.fd_rate), 
                0, minimum(oneweb_hist.fd_rate)] .* (sigtype.sig_freq/sig_freq)
    fd_rates_gps = [0, minimum(gps_hist.fd_rate)] .* (sigtype.sig_freq/sig_freq)
    fd_rates_iridium = [0, minimum(iridium_hist.fd_rate)] .* (sigtype.sig_freq/sig_freq)
    fd_rates_starlink = [0, minimum(starlink_hist.fd_rate)] .* (sigtype.sig_freq/sig_freq)
    fd_rates_oneweb = [0, minimum(oneweb_hist.fd_rate)] .* (sigtype.sig_freq/sig_freq)
    Pds = benchmark_acquisition(sigtype, Ts, CN0s, fd_rates, f_s; 
                                iterations=iterations, phi=phi, 
                                code_start_idx=code_start_idx,
                                f_if=f_if,
                                fd_ranges=fd_ranges)
    Pds_gps = Pds[:,:,1:2]
    Pds_iridium = Pds[:,:,3:4]
    Pds_starlink = Pds[:,:,5:6]
    Pds_oneweb = Pds[:,:,7:8]
    file = string(directory, file_names[i])
    save(file, 
         "sig_freq", sig_freq, 
         "Ts", Ts, 
         "CN0s", CN0s, 
         "fd_rates_gps", fd_rates_gps, 
         "fd_rates_iridium", fd_rates_iridium, 
         "fd_rates_starlink", fd_rates_starlink, 
         "fd_rates_oneweb", fd_rates_oneweb, 
         "Pds_gps", Pds_gps, 
         "Pds_iridium", Pds_iridium, 
         "Pds_starlink", Pds_starlink,
         "Pds_oneweb", Pds_oneweb, 
         "f_s", f_s, 
         "iterations", iterations)
end

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_names = ["ch4_section3_acquisition_l1ca_constellation_comparison",
              "ch4_section3_acquisition_l5q_constellation_comparison"]
items = [1, 2]
for i in items
    file_name = file_names[i]
    file = string(directory, file_name, ".jld")
    results = load(file, "Ts", "CN0s", "fd_rates_gps", "fd_rates_iridium",
                   "fd_rates_starlink", "fd_rates_oneweb", "Pds_gps",
                   "Pds_iridium", "Pds_starlink", "Pds_oneweb")
    Ts, CN0s, fd_rates_gps, fd_rates_iridium, fd_rates_starlink, 
    fd_rates_oneweb, Pds_gps, Pds_iridium, Pds_starlink, Pds_oneweb = results
    fig = figure(figsize=(8, 4.5))
    ax1 = plot_data_for_T(1, fig, Ts, CN0s, fd_rates_gps, Pds_gps,
                          plot_title="(a) GPS")
    legend()
    ax2 = plot_data_for_T(2, fig, Ts, CN0s, fd_rates_iridium, Pds_iridium;
                          plot_title="(b) Iridium")
    legend()
    ax3 = plot_data_for_T(3, fig, Ts, CN0s, fd_rates_starlink, Pds_starlink;
                          plot_title="(c) Starlink")
    ax4 = plot_data_for_T(4, fig, Ts, CN0s, fd_rates_oneweb, Pds_oneweb;
                          plot_title="(d) OneWeb")
    subplots_adjust(wspace=0.3, hspace=0.525, 
                    bottom=0.1, left=0.1, right=0.9, top=0.93)
    savefig(string(directory, "figures/", file_name, ".svg"), dpi=300)
end


