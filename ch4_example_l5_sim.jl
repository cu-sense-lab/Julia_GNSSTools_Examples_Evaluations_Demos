using GNSSTools
using Statistics
using Random
using StatsBase
using LinearAlgebra
using Base.Threads
using JLD
using PyPlot
pygui(true)
include("ch4_data_plot_functions.jl")
include("ch4_constellation_hist_structs.jl")
include("ch4_simulation_function.jl")


directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
folder = "figures/l5_sim/"

# doppler_file_name = string(directory, "ch1_doppler_curves.jld")
# doppler_raw_gps, doppler_rate_raw_gps, sig_freq = load(doppler_file_name, 
#                                                        "doppler_raw_gps",
#                                                        "doppler_rate_raw_gps",
#                                                        "freq")
# gps_hist = initiate_hist(doppler_raw_gps, doppler_rate_raw_gps; nbins=1000)
#
# t_length = 2
# channel = "Q"
# f_s = 25e6
# f_if = 0
# Tsys = 535
# nADC = 8
# N = floor(Int, f_s*20e-3)
# B = 25e6
# prn_num = 9
# upsample_factor = 4
# prns = []
# for i in 1:prn_num
#     prn = rand(1:32)
#     while prn in prns
#         prn = rand(1:32)
#     end
#     push!(prns, prn)
# end
# prns = sort(prns)
# phis = rand(0:0.0001:2π, prn_num)
# cn0s = rand(35:0.01:45, prn_num)
# f_ds = Array{Float64}(undef, prn_num)
# fd_rates = Array{Float64}(undef, prn_num)
# for i in 1:prn_num
#     f_d, fd_rate = sample_distribution(gps_hist)
#     f_ds[i] = f_d
#     fd_rates[i] = fd_rate
# end
# code_start_idxs = rand(1:N, prn_num) + Rational.(1, (rand(1:upsample_factor, prn_num)))
# h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # Rakon IT5300B TCXO

# save(string(directory, folder, "ch4_l5_sim_parms.jld"),
#             "t_length", t_length,
#             "f_s", f_s,
#             "f_if", f_if,
#             "Tsys", Tsys,
#             "nADC", nADC,
#             "N", N,
#             "B", B,
#             "prn_num", prn_num,
#             "upsample_factor", upsample_factor,
#             "prns", prns,
#             "phis", phis,
#             "cn0s", cn0s,
#             "f_ds", f_ds,
#             "fd_rates", fd_rates,
#             "code_start_idxs", code_start_idxs,
#             "h_parms", h_parms,
#             "channel", channel)

parms = load(string(directory, folder, "ch4_l5_sim_parms.jld"),
             "t_length",
             "f_s",
             "f_if",
             "Tsys",
             "nADC",
             "N",
             "B",
             "prn_num",
             "upsample_factor",
             "prns",
             "phis",
             "cn0s",
             "f_ds",
             "fd_rates",
             "code_start_idxs",
             "h_parms",
             "channel")

t_length, f_s, f_if, Tsys, nADC, N, B, prn_num, upsample_factor, prns,
phis, cn0s, f_ds, fd_rates, code_start_idxs, h_parms, channel = parms

sigtype = define_l5_code_type(;channel="both", B=B)

data = simulate_signal(sigtype, f_s, t_length, prns, cn0s.+3, f_ds, fd_rates, 
                       phis, code_start_idxs, h_parms; 
                       nADC=nADC, f_if=f_if, Tsys=Tsys)

sigtype = define_l5_code_type(;channel=channel, B=B)
replica = definereplica(sigtype, 25e6, 20e-3)
detectable_prns = []
for i in 1:32
    fd_est, n0_est, SNR_est, 
    P_d, above_threshold = courseacquisition(data, replica, i; M=1, 
                                             Δfd=1/(1*replica.t_length))
    if above_threshold && (P_d >= 0.9)
        detectable = true
        push!(detectable_prns, i)
    else
        detectable = false
    end
    f_code_d, f_code_dd = GNSSTools.calc_doppler_code_rate(nh_chipping_rate, 
                                                           replica.signal_type.sig_freq, 
                                                           fd_est, 0)
    n0 = calcinitcodephase(length(nh20), f_code_d, f_code_dd, data.f_s, 
                           n0_est)
    println("PRN $i\t$SNR_est\t$fd_est\t$n0_est\t$n0\t$(Float64(round(P_d, digits=2)))\t$detectable")
end
println("$(length(detectable_prns)) Detected PRNs: $detectable_prns")

fig = figure(figsize=(7.5, 5))
acquisition_T = 20e-3
fine_acq_T = 100e-3
fine_acq_sub_T = 1e-3
tracking_T = 20e-3
M = 1
linewidth = 0.5
tracking_σᵩ = Array{Float64}(undef, length(detectable_prns))
tracking_snr = Array{Float64}(undef, length(detectable_prns))
tracking_code_phase_err = Array{Float64}(undef, length(detectable_prns))
for i in 1:length(detectable_prns)
    prn = detectable_prns[i]
    ax = fig.add_subplot(3, 3, i)
    tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    title("PRN $prn")
    acqresults, trackresults = process(data, sigtype, prn, "Q";
                                       show_plot=false,
                                       fine_acq_method=:fft, 
                                       state_num=3,
                                       dynamickf=true, 
                                    #    h₀=h_parms[3], h₋₂=h_parms[1], 
                                       dll_b=1,
                                       acquisition_T=acquisition_T, 
                                       fine_acq_T=fine_acq_T,
                                       fine_acq_sub_T=fine_acq_sub_T,
                                       tracking_T=tracking_T,
                                       use_fine_acq=true,
                                       M=M)
    ax.plot(trackresults.t, real.(trackresults.ZP), linewidth=linewidth)
    ax.plot(trackresults.t, imag.(trackresults.ZP), linewidth=linewidth)
    N = floor(Int, length(trackresults.t)/2)
    tracking_σᵩ[i] = std(trackresults.dphi_meas[N:end]) * 180 / π
    tracking_snr[i] = mean(trackresults.SNR[N:end])
    tracking_code_phase_err[i] = std(trackresults.code_err_meas[N:end])
end
subplots_adjust(hspace=0.4, wspace=0.4,  top=0.88, left=0.02, right=0.98, 
                bottom=0.02)
plot_title = string("Course Acq. T=$(floor(Int, acquisition_T*1000))ms    ",
                    "Fine Acq. T=$(floor(Int, fine_acq_T*1000))ms    ",
                    "Tracking T=$(floor(Int, tracking_T*1000))ms    ",
                    "M=$(M)")
suptitle(plot_title)
T = floor(Int, tracking_T*1e3)
savefig(string(directory, folder, "ch4_l5_T=$(T)ms_acquired_prns.svg"), dpi=300)


T = 20e-3
for i in 1:length(detectable_prns)
    detected_prn = detectable_prns[i]
    j = argmax(detected_prn .== prns)
    acqresults, trackresults, corr_result, SNR_est, 
    P_d, above_threshold = process(data, 
                                sigtype, 
                                detected_prn, channel;
                                fine_acq_method=:fft, 
                                return_corrresult=true,
                                return_Pd=true,
                                state_num=3, 
                                #    h₀=h_parms[3], 
                                #    h₋₂=h_parms[1],
                                q_a=1,
                                σω=1,
                                dll_b=1,
                                acquisition_T=20e-3, 
                                fine_acq_T=100e-3,
                                fine_acq_sub_T=1e-3,
                                tracking_T=T,
                                show_plot=false,
                                M=1,
                                use_fine_acq=true);
    # Calculate truth info
    dt = T
    code_start_idx = code_start_idxs[j]
    f_d = f_ds[j]
    fd_rate = fd_rates[j]
    phi = phis[j]
    cn0 = cn0s[j]
    t, doppler, doppler_rate, phi, code_phase, snr = get_signal_truth_parms(t_length, 
                                                                        f_s, 
                                                                        code_start_idx, 
                                                                        f_d, fd_rate, 
                                                                        phi, cn0, Tsys, 
                                                                        sigtype, channel, 
                                                                        T; dt=dt, B=B)
    name = "pdfs/ch4_l5$(lowercase(channel))_T=$(floor(Int, T*1e3))ms_prn_$(trackresults.prn).pdf"
    file_name = string(directory, folder, name)
    plot_example_track_results(trackresults; f_if=0, code_length=l1ca_code_length,
                            save_to_file=file_name,
                            bottom=0.075, top=0.95, left=0.125, right=0.875, 
                            hspace=0.75, wspace=0.35, marker_size=1,
                            line_width=0.8,
                            truth_t=t, truth_doppler=doppler, 
                            truth_doppler_rate=doppler_rate,
                            truth_code_phase=code_phase,
                            truth_snr=snr, truth_phi=phi)
end

