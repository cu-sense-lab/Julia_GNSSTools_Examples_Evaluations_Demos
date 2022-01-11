using GNSSTools
using Statistics
using ProgressMeter
using JLD
using PyPlot
pygui(true)

h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
CN0 = 45.  # dB⋅Hz
f_d = 0
nADC = 16
code_start_idx = 1000
iterations = 100
Ts = Array(1e-3:1e-3:100e-3)
# t_length = maximum(Ts)

# L1 C/A
f_s = 2.048e6
sigtype = define_l1ca_code_type(; sig_freq=1*L1_freq)
# generatesignal!(signal)
l1ca_SNRs = zeros(length(Ts))
p = Progress(length(Ts)*iterations, 1, "Processing L1 C/A...")
for i in 1:length(Ts)
    SNR = 0
    T = Ts[i]
    signal = definesignal(sigtype, f_s, T; 
                          prn=1, code_start_idx=code_start_idx, f_d=f_d, fd_rate=0, 
                          CN0=CN0, include_phase_noise=false, f_if=1.25e6,
                          include_thermal_noise=true, include_adc=false, nADC=nADC,
                          include_carrier=true, receiver_h_parms=h_parms, 
                          skip_noise_generation=false)
    replica = definesignal(sigtype, f_s, T; f_d=f_d, allocate_noise_vectors=false, 
                           skip_noise_generation=true)
    generatereplica!(replica)
    for j in 1:iterations
        # definesignal!(signal; new_thermal_noise=true, new_phase_noise=false)
        # generatesignal!(signal)
        generatereplica!(signal)
        A = sqrt(2*GNSSTools.k*535*10^(CN0/10))
        σₙ = sqrt(GNSSTools.k*2.046e6*535)
        signal.data .= A.*signal.data .+ σₙ.*randn(signal.sample_num)
        corr_result = fft_correlate(signal.data[1:replica.sample_num], replica.data)
        corr_result_fft = fft(corr_result)
        PS = 2*maximum(abs2.(corr_result_fft))
        # PS = abs2(corr_result[code_start_idx])
        PN = (sum(abs2.(corr_result_fft)) - PS) / (length(corr_result_fft) - 2)
        SNR += PS/PN
        next!(p)
    end
    l1ca_SNRs[i] = 10log10(SNR/iterations)
end


fig = figure()
plot(Ts.*1000, l1ca_SNRs, label="L1 C/A")
xlabel("Integration Times (ms)")
ylabel("Acquired SNR (dB)")
legend()
xscale("log")

save("ch2_snr_vs_integration_time.jld", "l1ca_SNRs", l1ca_SNRs)
savefig("figures/ch2_snr_vs_integration_time.pdf", dpi=300)

