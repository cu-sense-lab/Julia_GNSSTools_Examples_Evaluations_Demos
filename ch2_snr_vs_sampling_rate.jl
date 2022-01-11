using GNSSTools
using Statistics
using ProgressMeter
using PyPlot
pygui(true)

h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
CN0 = 45.  # dB⋅Hz
f_d = 0
iterations = 100
T = 1e-3
f_s₀ = 2.048e6
f_ss = Array(1:100)
t_length = T

# L1 C/A
sigtype = define_l1ca_code_type(; sig_freq=1*L1_freq)
# generatesignal!(signal)
l1ca_SNRs = zeros(length(f_ss))
p = Progress(length(f_ss)*iterations, 1, "Processing L1 C/A...")
for i in 1:length(f_ss)
    SNR = 0
    f_s = f_ss[i] * f_s₀
    replica = definesignal(sigtype, f_s, T; f_d=f_d, allocate_noise_vectors=false, 
                           skip_noise_generation=true)
    generatereplica!(replica)
    signal = definesignal(sigtype, f_s, t_length; 
                          prn=1, code_start_idx=1000, f_d=f_d, fd_rate=0, 
                          CN0=CN0, include_phase_noise=false, 
                          include_thermal_noise=true, include_adc=false, 
                          include_carrier=true, receiver_h_parms=h_parms,
                          skip_noise_generation=true)
    for j in 1:iterations
        definesignal!(signal; new_thermal_noise=true, new_phase_noise=false)
        generatesignal!(signal)
        corr_result = fft_correlate(signal.data[1:replica.sample_num], replica.data)
        PS = maximum(abs2.(corr_result))
        PN = ((sum(abs2.(corr_result))) - PS) / (replica.sample_num - 1)
        SNR += PS/PN
        next!(p)
    end
    l1ca_SNRs[i] = 10log10(SNR/iterations)
end


fig = figure()
plot(f_ss.*f_s₀./1e6, l1ca_SNRs, label="L1 C/A")
# plot(Ts, l5_SNRs, label="L5")
xlabel("Sampling Rate (MHz)")
ylabel("Acquired SNR (dB)")
legend()
xscale("log")
savefig("figures/ch2_snr_vs_sampling_rate.pdf", dpi=300)
