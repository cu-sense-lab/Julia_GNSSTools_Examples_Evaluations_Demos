using GNSSTools
using Statistics
using Random
using StatsBase
using LinearAlgebra
using ProgressMeter
using JLD
using PyPlot
pygui(true)

struct ConstellationHistogram{T1,T2}
    nbins::Int
    dopplers::Array{Float64,1}
    doppler_rates::Array{Float64,1}
    histogram::T1
    weights::T2
    fd_rate::Vector{Float64}
    f_d::Vector{Float64}
end


function initiate_hist(dopplers, doppler_rates; nbins=1000)
    histogram = normalize(fit(Histogram, (doppler_rates, dopplers), 
                              nbins=nbins), mode=:probability)
    weights = collect(Iterators.flatten(histogram.weights))
    weights = Weights(weights)
    fd_rate, f_d = meshgrid(Array(histogram.edges[1][1:end-1]), 
                            Array(histogram.edges[2][1:end-1]))
    fd_rate = collect(Iterators.flatten(fd_rate))
    f_d = collect(Iterators.flatten(f_d))
    return ConstellationHistogram(nbins, dopplers, doppler_rates, 
                                  histogram, weights, fd_rate, f_d)
end


function sample_distribution(histogram::ConstellationHistogram) 
    f_d = sample(histogram.f_d, histogram.weights)
    fd_rate = sample(histogram.fd_rate, histogram.weights)
    return (f_d, fd_rate)
end


# "doppler_gps",
# "doppler_iridium",
# "doppler_starlink",
# "doppler_oneweb",
# "doppler_rate_gps",
# "doppler_rate_iridium",
# "doppler_rate_starlink",
# "doppler_rate_oneweb",
# "doppler_bounds_gps",
# "doppler_bounds_iridium",
# "doppler_bounds_starlink",
# "doppler_bounds_oneweb",
# "doppler_raw_gps",
# "doppler_raw_iridium",
# "doppler_raw_starlink",
# "doppler_raw_oneweb",
# "doppler_rate_raw_gps",
# "doppler_rate_raw_iridium",
# "doppler_rate_raw_starlink",
# "doppler_rate_raw_oneweb",
# "freq",
# "user_lla",
# "t_range",
# "min_elevation")


function benchmark_tracking(sigtype::GNSSTools.SignalType, 
                            histogram::ConstellationHistogram, 
                            Ts, CN0s, fd_rates, f_s; iterations=100,
                            include_phase_noise=true, 
                            include_thermal_noise=true,
                            include_adc=true, 
                            include_carrier=true, 
                            h_parms=[1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2,
                            code_start_idx=missing, phi=missing,
                            nADC=8, fd_range=5000, f_if=1.25e6, prn=1, M=20,
                            fd_center=0, Pd_threshold=0.9, p_fa=1e-7,
                            acquisition_T=1e-3, fine_acq_T=10e-3,
                            fine_acq_method=:fft, use_fine_acq=true,
                            h₀=h_pamrs[3], h₋₂=h_parms[1], dll_b=0.1,
                            state_num=2, q_a=1, t_length=1, q_mult=1,
                            cov_mult=1, R_mult=1, dynamickf=true) 
    if sigtype.include_I && sigtype.include_Q
        B = max(sigtype.B_I, sigtype.B_Q)
    elseif sigtype.include_I && !sigtype.include_Q
        B = sigtype.B_I
    elseif !sigtype.include_I && sigtype.include_Q
        B = sigtype.B_Q
    else
        error("No CodeType defined in SigType struct.")
    end
    signal = definesignal(sigtype, f_s, t_length; 
                          prn=prn, f_if=f_if, nADC=nADC,
                          include_phase_noise=include_phase_noise, 
                          include_thermal_noise=include_thermal_noise, 
                          include_adc=include_adc, 
                          include_carrier=include_carrier, 
                          receiver_h_parms=h_parms)
    Pds = Array{Float64}(undef, length(Ts), length(CN0s))
    p = Progress(length(Ts)*length(CN0s)*iterations, 1, "Processing...")
    for i in 1:length(Ts)
        tracking_T = Ts[i]
        N_over_2 = floor(Int, t_length/tracking_T/2)
        for j in 1:length(CN0s)
            CN0 = CN0s[j]
            
            for iteration in 1:iterations
                if ismissing(phi)
                    ϕ₀_truth = rand(0:0.0001:2π)
                else
                    ϕ₀_truth = phi
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
                f_d, fd_rate = sample_distribution(histogram)
                definesignal!(signal; CN0=CN0, phi=ϕ₀_truth, f_d=f_d,
                              code_start_idx=n0_truth, fd_rate=fd_rate,
                              new_phase_noise=new_phase_noise, 
                              new_thermal_noise=new_thermal_noise)
                generatesignal!(signal)
                acqresults, trackresults, corr_result, SNR_est, 
                P_d, above_threshold = process(signal, 
                                               sigtype, 
                                               prn, "I"; fd_range=fd_range, 
                                               fine_acq_method=fine_acq_method, 
                                               return_corrresult=true, 
                                               return_Pd=true,
                                               h₀=h₀, 
                                               h₋₂=h₋₂, 
                                               q_mult=q_mult,
                                               cov_mult=cov_mult,
                                               R_mult=R_mult,
                                               q_a=q_a,
                                               dynamickf=dynamickf,
                                               state_num=state_num,
                                               dll_b=dll_b,
                                               acquisition_T=acquisition_T,
                                               fine_acq_T=fine_acq_T,
                                               tracking_T=tracking_T,
                                               show_plot=false,
                                               M=M);
                if P_d > Pd_threshold
                    dphi_meas_std = std(trackresults.dphi_meas[N_over_2:end])
                    code_err_meas_std = std(trackresults.code_err_meas[N_over_2:end])
                    P = mean.(trackresults.P[N_over_2:end])
                end
                next!(p)
            end
        end
    end
end


file_name = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/ch1_doppler_curves.jld"
doppler_raw_gps, doppler_rate_raw_gps, sig_freq = load(file_name, 
                                                       "doppler_raw_gps",
                                                       "doppler_rate_raw_gps",
                                                       "freq")


gps_hist = initiate_hist(doppler_raw_gps, doppler_rate_raw_gps; nbins=1000)
# f_d, fd_rate = sample_distribution(gps_hist)
f_s = 2.5e6  # Hz
Ts = [1, 2, 5, 10] .* 1e-3  # seconds
CN0s = [25, 30, 35, 40, 45, 50]  # dB
fd_rates_gps = unique(gps_hist.fd_rate)[1:30:end]  # Hz/second
sigtype = define_l1ca_code_type(maximum(Ts))


