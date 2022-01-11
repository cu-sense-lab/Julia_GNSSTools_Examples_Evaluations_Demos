using GNSSTools
using Statistics
using ProgressMeter
using PyPlot
pygui(true)

# h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
t_length = 1e-3
CN0s = Array(30:0.25:50)
iterations = 1000
T = 1e-3
f_d = 0

# L1 C/A
f_s = 2.048e6
sigtype = define_l1ca_code_type(; sig_freq=1*L1_freq)
replica = definesignal(sigtype, f_s, T)
signal = definesignal(sigtype, f_s, t_length; 
                      prn=1, code_start_idx=1000, f_d=f_d, fd_rate=0, 
                      CN0=36, include_phase_noise=false, 
                      include_thermal_noise=true, include_adc=true, 
                      include_carrier=true, receiver_h_parms=h_parms)
l1ca_SNRs = zeros(length(CN0s))
p = Progress(length(CN0s)*iterations, 1, "Processing L1 C/A...")
for i in 1:length(CN0s)
    SNR = 0
    for j in 1:iterations
        definesignal!(signal; CN0=CN0s[i], new_thermal_noise=true, 
                      new_phase_noise=false)
        generatesignal!(signal)
        fd_est, n0_est, SNR_est = courseacquisition(signal, replica, 1)
        SNR += (10^(SNR_est/10))^2
        next!(p)
    end
    l1ca_SNRs[i] = 10log10(SNR/iterations)
end

# save("ch4_acquisition_performance.jld", "results", results)


# L%I
f_s = 10*2.048e6
sigtype = definesignaltype(definecodetype(l5i_codes, L5_chipping_rate), 
                           L5_freq, "I")
replica = definesignal(sigtype, f_s, T)
signal = definesignal(sigtype, f_s, t_length; 
                      prn=1, code_start_idx=1000, f_d=f_d, fd_rate=0,  
                      CN0=36, include_phase_noise=false, 
                      include_thermal_noise=true, include_adc=true, 
                      include_carrier=true, receiver_h_parms=h_parms)
l5_SNRs = zeros(length(CN0s))
p = Progress(length(CN0s)*iterations, 1, "Processing L5...")
for i in 1:length(CN0s)
    SNR = 0
    for j in 1:iterations
        definesignal!(signal; CN0=CN0s[i], new_thermal_noise=true, 
                      new_phase_noise=false)
        generatesignal!(signal)
        fd_est, n0_est, SNR_est = courseacquisition(signal, replica, 1)
        SNR += (10^(SNR_est/10))^2
        next!(p)
    end
    l5_SNRs[i] = 10log10(SNR/iterations)
end


fig = figure()
plot(CN0s, l1ca_SNRs, label="L1 C/A")
plot(CN0s, l5_SNRs, label="L5")
xlabel("C/N₀ (dB⋅Hz)")
ylabel("Acquired SNR (dB)")
legend()
savefig("figures/ch2_l1ca_vs_l5_snr_at_2x_sampling_rate.pdf", dpi=300)