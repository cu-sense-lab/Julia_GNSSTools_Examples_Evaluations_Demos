using GNSSTools
using Statistics
using ProgressMeter
using PyPlot
pygui(true)


CN0s = [25, 30, 35, 40]  # dB⋅Hz
h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
# h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16]./10e6^2   # OCXO

t_length = 5.  # seconds
wait_time = 1.  # seconds (the time to wait before calculating STD of phase)
f_s = 5e6  # Hz
sig_freq = L1_freq  # Hz
f_d = 0.  # Hz
f_if = 0.  # Hz
fd_rates = 0:-5:-350  # Hz/s
integration_times = range(1e-3, 40e-3, step=1e-3)
iterations = 1
acquisition_T = 1e-3  # seconds
fine_acq_T = 10e-3  # seconds
fine_acq_method = :fft
prn = 1

sigtype = define_l1ca_code_type(;sig_freq=sig_freq)
signal = definesignal(sigtype, f_s, t_length; prn=prn, f_d=f_d, f_if=f_if)

results = Array{Float64}(undef, length(CN0s), length(fd_rates), length(integration_times))

p = Progress(iterations*length(CN0s)*length(fd_rates)*length(integration_times), 1, "Processing...")
for i in 1:length(CN0s)
    for j in 1:length(fd_rates)
        for k in 1:length(integration_times)
            Δϕσ = 0.
            for iteration in 1:iterations
                CN0 = CN0s[i]
                fd_rate = fd_rates[j]
                T = integration_times[k]
                ϕ₀ = rand(0:0.0001:2π)
                code_start_idx = rand(1:signal.sample_num)
                definesignal!(signal; fd_rate=fd_rate, CN0=CN0,
                              code_start_idx=code_start_idx, phi=ϕ₀,
                              new_thermal_noise=true, new_phase_noise=true)
                generatesignal!(signal)
                # println([i,j,k])
                acqresults, trackresults = process(signal, sigtype, prn,
                                                   acquisition_T=acquisition_T,
                                                   fine_acq_T=fine_acq_T,
                                                   fine_acq_method=fine_acq_method,
                                                   tracking_T=T, show_plot=false,
                                                   h₋₂=h_parms[1], h₀=h_parms[3])
                dphi_start_idx = floor(Int, wait_time/T)
                Δϕσ += std(trackresults.dphi_meas[dphi_start_idx:end])*180/π
                next!(p)
            end
            results[i,j,k] = Δϕσ/iterations
        end
    end
end
