using GNSSTools
using Statistics
using Base.Threads
using ProgressMeter
using JLD
using FFTW
using PyPlot
pygui(true)
include("ch4_constellation_hist_structs.jl")

# CN0s = [30, 35, 40, 45]  # dB⋅Hz
# h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
# h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16]./10e6^2   # OCXO

# f_s = 5e6  # Hz
# sig_freq = L1_freq  # Hz
# f_d = 0.  # Hz
# f_if = 0.  # Hz
# fd_rates = 0:-10:-350  # Hz/s
# integration_times = range(1e-3, 20e-3, step=1e-3)
# t_length = maximum(integration_times)  # seconds
# iterations = 1
# prn = 1

# sigtype = define_l1ca_code_type(;sig_freq=sig_freq)
# signal = definesignal(sigtype, f_s, t_length; prn=prn, f_d=f_d, f_if=f_if,
#                       include_phase_noise=false)

# results = Array{Float64}(undef, length(CN0s), length(fd_rates), length(integration_times))

# p = Progress(iterations*length(CN0s)*length(fd_rates)*length(integration_times), 1, "Processing...")
# for i in 1:length(CN0s)
#     for j in 1:length(fd_rates)
#         for k in 1:length(integration_times)
#             SNR = 0.
#             for iteration in 1:iterations
#                 CN0 = CN0s[i]
#                 fd_rate = fd_rates[j]
#                 T = integration_times[k]
#                 ϕ₀ = rand(0:0.0001:2π)
#                 code_start_idx = rand(1:signal.sample_num)
#                 definesignal!(signal; fd_rate=fd_rate, CN0=CN0,
#                               code_start_idx=code_start_idx, phi=ϕ₀,
#                               new_thermal_noise=true, new_phase_noise=false)
#                 generatesignal!(signal)
#                 # println([i,j,k])
#                 replica = definesignal(sigtype, f_s, T;
#                                        skip_noise_generation=true,
#                                        allocate_noise_vectors=false)
#                 fd_course, n0_est, SNR_est = courseacquisition(signal, 
#                                                                replica,
#                                                                prn)
#                 SNR += 10^(SNR_est/10)
#                 next!(p)
#             end
#             results[i,j,k] = 10*log10(SNR/iterations)
#         end
#     end
# end

# save("ch4_acquisition_performance_l1ca.jld", "results", results)

# results = load("ch4_acquisition_performance.jld", "results")

# max_val = maximum(results)
# min_val = minimum(results)
# heatmap_extent = [1e3*integration_times[1], 1e3*integration_times[end], fd_rates[1], fd_rates[end]]
# fig = figure(figsize=(7.5,2))
# for i in 1:length(CN0s)
#     ax = fig.add_subplot(1, length(CN0s), i)
#     global h = ax.imshow(results[i,:,:], aspect="auto", vmin=min_val, vmax=max_val, extent=heatmap_extent)
#     xlabel("Integration Time (ms)")
#     ylabel("Doppler Rate (Hz/s)")
#     title("C/N₀: $(CN0s[i])dB⋅Hz")
# end
# cbar_ax = fig.add_axes([0.92, 0.15, 0.04, 0.7])
# cbar = fig.colorbar(h, cax=cbar_ax)
# cbar.set_label("SNR (dB)")
# subplots_adjust(hspace=0.3, wspace=0.4, top=0.91, left=0.05, right=0.9, bottom=0.12)
# savefig("figures/ch4_acquisition_performance.pdf", dpi=300)


################################################################################
CN0s = [45]  # dB⋅Hz
h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
# h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16]./10e6^2   # OCXO


directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"

integration_times = [1e-3, 2e-3, 5e-3, 10e-3]  # seconds
sigtype = define_l1ca_code_type(maximum(integration_times))
f_s = 2.5e6  # Hz
sig_freq = L1_freq  # Hz
file_name = "ch4_acquisition_performance_nonheatmap_l1ca"


# sigtype = define_l5_code_type(;channel="Q")
# integration_times = [1e-3, 20e-3, 40e-3, 100e-3]  # seconds
# f_s = 25e6  # Hz
# sig_freq = L5_freq  # Hz
# file_name = "ch4_acquisition_performance_nonheatmap_l5"

f_d = 0.  # Hz 
f_if = 0.  # Hz
# fd_rate_max = 10e3 * (sig_freq/L1_freq)
fd_rate_max = 10e3 * (sig_freq/L1_freq)
fd_rates = 0:-10:-fd_rate_max  # Hz/s
iterations = 100
prn = 1

results = Array{Float64}(undef, length(CN0s), length(fd_rates), length(integration_times))

p = Progress(iterations*length(CN0s)*length(fd_rates)*length(integration_times), 1, "Processing...")
@threads for k in 1:length(integration_times)
    T = integration_times[k]
    signal = definesignal(sigtype, f_s, T; prn=prn, f_if=f_if,
                      include_phase_noise=true, include_thermal_noise=true,
                      nADC=8,
                      receiver_h_parms=h_parms)
    replica = definereplica(sigtype, f_s, T; f_d=0)
    generatereplica!(replica)
    for i in 1:length(CN0s)
        for j in 1:length(fd_rates)
                SNR = 0.
                CN0 = CN0s[i]
                fd_rate = fd_rates[j]
                for iteration in 1:iterations
                    ϕ₀ = rand(0:0.0001:2π)
                    # f_d = rand(-500:0.1:500)
                    f_d = 0
                    code_start_idx = rand(1:signal.sample_num)
                    definesignal!(signal; fd_rate=fd_rate, CN0=CN0,
                                code_start_idx=code_start_idx, phi=ϕ₀,
                                new_thermal_noise=true, new_phase_noise=true,
                                f_d=f_d)
                    generatesignal!(signal)
                    corr_result = fft_correlate(signal.data[1:replica.sample_num], 
                                                replica.data)
                    corr_result = abs2.(corr_result)
                    PS = maximum(corr_result)
                    PN = (sum(corr_result) - PS) / (length(corr_result) - 1)
                    SNR += PS/(2*PN)
                    next!(p)
                end
                results[i,j,k] = 10log10(SNR/iterations)
        end
    end
end

save(string(directory, file_name, ".jld"), "results", results,
     "integration_times", integration_times, "fd_rates", fd_rates,
     "sig_freq", sig_freq, "f_d", f_d, "f_if", f_if, "f_s", f_s)

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


const_names = ["GPS", "Iridium", "Starlink", "OneWeb"]
const_max_fd_rates = [minimum(doppler_rate_raw_gps),
                      minimum(doppler_rate_raw_iridium),
                      minimum(doppler_rate_raw_starlink),
                      minimum(doppler_rate_raw_oneweb)] 

file_names = ["ch4_acquisition_performance_nonheatmap_l1ca",
              "ch4_acquisition_performance_nonheatmap_l5"]
fig = figure(figsize=(7.5, 3.5))
for j in 1:length(file_names)
    file_name = file_names[j]
    results, fd_rates, integration_times = load(string(directory, file_name, ".jld"), 
                                                "results",
                                                "fd_rates",
                                                "integration_times")
    ax = fig.add_subplot(1,2,j)
    for i in 1:length(integration_times)
        plot(abs.(fd_rates), results[1,:,i], 
            label="T = $(Int(integration_times[i]*1e3))ms")
    end
    for k in 1:length(const_names)
        axvline(x=abs(const_max_fd_rates[k]), c="grey", linestyle=":")
    end
    xlabel("Doppler Rate (Hz/s)")
    ylabel("SNR (dB)")
    xscale("log")
    legend()
    letter = ('a' + (j-1))
    title("($(letter))")
end
subplots_adjust(wspace=0.3, bottom=0.15, left=0.08, right=0.92, top=0.93)
savefig(string(directory, "figures/ch4_acquisition_high_doppler.svg"), dpi=300)
