using GNSSTools
using JLD
using PyPlot
pygui(true)

function vsm(trackresults, K)
    ZP = trackresults.ZP
    M = floor(Int, length(ZP)/K)
    CN0_est = Array{Float64}(undef, M)
    t_est = Array{Float64}(undef, M)
    for i in 1:M
        t_est[i] = trackresults.t[i*K]
        zp = ZP[(i-1)*K+1:i*K]
        Zₖ = real.(zp).^2 + imag.(zp).^2
        Z̄ = mean(Zₖ)
        σ² = var(Zₖ)
        A = sqrt(Z̄^2 - σ²)
        σ²_iq = (Z̄ - A)/2
        CN0 = 10log10(A/(2*σ²_iq*trackresults.T))
        CN0_est[i] = CN0
    end
    return (CN0_est, t_est)
end

function snr_est(trackresults, K)
    ZP = trackresults.ZP
    M = floor(Int, length(ZP)/K)
    SNR_est = Array{Float64}(undef, M)
    t_est = Array{Float64}(undef, M)
    for i in 1:M
        t_est[i] = trackresults.t[i*K]
        zp = ZP[(i-1)*K+1:i*K]
        I² = mean(real.(zp))^2
        Q² = mean(imag.(zp))^2
        SNR_est[i] = I² / (2*Q²)
    end
    return (SNR_est, t_est)
end


h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # Rakon IT5300B TCXO


#------------------- Lab Data Sample for L1 C/A Signal ------------------------#
# file = "/home/sjbilardi/Documents/SeNSeLab/GNSSTools.jl/test/data/known_oscillator_data/20210327_013800_000005_Aero-North-TCXO_1575.42_25.sc16"
# file = "/home/sjbilardi/Documents/SeNSeLab/GNSSTools.jl/test/data/20201217_224100_GPSL1_1575.42_25.0.sc8"
file = "/home/sjbilardi/Documents/SeNSeLab/GNSSTools.jl/test/data/20201211_211745_GPSL1_1575.42_25.0.sc4"
data = loaddata(:sc4, file, 25e6, L1_freq, L1_freq, 3, skip_to=1)

replica = definereplica(define_l1ca_code_type(), 25e6, 2e-3)
detectable_prns = []
for i in 1:32
    fd_est, n0_est, SNR_est, 
    P_d, above_threshold = courseacquisition(data, replica, i; M=20, 
                                             Δfd=1/(1*replica.t_length))
    if (P_d >= 0.9)
        detectable = true
        push!(detectable_prns, i)
    else
        detectable = false
    end
    println("PRN $i\t$SNR_est\t$(Float64(round(P_d, digits=2)))\t$detectable")
end
println("$(length(detectable_prns)) Detected PRNs: $detectable_prns")

fig = figure()
acquisition_T = 2e-3
fine_acq_T = 10e-3
tracking_T = 1e-3
M = 10
for i in 1:32
    ax = fig.add_subplot(8, 4, i)
    tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    title("PRN $i")
    acqresults, trackresults = process(data, define_l1ca_code_type(), i;
                                       show_plot=false,
                                       fine_acq_method=:fft, state_num=2,
                                       dynamickf=true, 
                                       # h₀=h_parms[3], 
                                       # h₋₂=h_parms[1], 
                                       dll_b=0.01,
                                       acquisition_T=acquisition_T, 
                                       fine_acq_T=fine_acq_T,
                                       tracking_T=tracking_T,
                                       M=M)
    ax.plot(trackresults.t, real.(trackresults.ZP))
    ax.plot(trackresults.t, imag.(trackresults.ZP))
end
subplots_adjust(hspace=0.4, wspace=0.4,  top=0.93, left=0.08, right=0.93, 
                bottom=0.01)
plot_title = string("Course Acq. T=$(floor(Int, acquisition_T*1000))ms;    ",
                    "Fine Acq. T=$(floor(Int, fine_acq_T*1000))ms;    ",
                    "Tracking T=$(floor(Int, tracking_T*1000))ms;    ",
                    "Noncoherent Integrations=$(M)")
suptitle(plot_title)


acqresults, trackresults, corr_result, SNR_est, 
P_d, above_threshold = process(data, 
                               define_l1ca_code_type(), 
                               3;
                               fine_acq_method=:fft, 
                               return_corrresult=true,
                               return_Pd=true,
                               state_num=3, 
                               # h₀=h_parms[3], 
                               # h₋₂=h_parms[1],
                               q_a=1,
                               dll_b=1,
                               acquisition_T=2e-3, 
                               fine_acq_T=10e-3,
                               tracking_T=1e-3,
                               M=10,
                               use_fine_acq=true,
                               show_plot=false,
                               σω=1);

#------------------------------------------------------------------------------#


#--------------------- Lab Data Sample for L5 Signal --------------------------#
# file = "/home/sjbilardi/Documents/SeNSeLab/GNSSTools.jl/test/data/20201211_211300_GPSL5_1176.45_25.0.sc4"
file = "/home/sjbilardi/Documents/SeNSeLab/GNSSTools.jl/test/data/20200622_203230_000100_g5_ref.jld"
# data = loaddata(:sc4, file, 25e6, L5_freq, L5_freq, 1, skip_to=0)
data = load(file, "data")

sigtype = definesignaltype(definecodetype([l5q_codes, nh20], 
                                          [L5_chipping_rate, nh_chipping_rate]),
                                          L5_freq, "Q")
replica = definereplica(sigtype, 25e6, 20e-3)
detectable_prns = []
for i in 1:32
    fd_est, n0_est, SNR_est, 
    P_d, above_threshold = courseacquisition(data, replica, i; M=1, 
                                             Δfd=1/(1*replica.t_length))
    if above_threshold && (P_d >= 0.9)
        push!(detectable_prns, i)
    end
    println("PRN $i\t$SNR_est\t$(Float64(round(P_d, digits=2)))\t$above_threshold")
end
println("Detected PRNs: $detectable_prns")

fig = figure()
acquisition_T = 1e-3
fine_acq_T = 100e-3
tracking_T = 1e-3
M = 20
for i in 1:32
    ax = fig.add_subplot(8, 4, i)
    tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    title("PRN $i")
    acqresults, trackresults = process(data, sigtype, i, "Q";
                                       show_plot=false,
                                       fine_acq_method=:fft, state_num=2,
                                       dynamickf=true, h₀=h_parms[3], 
                                       h₋₂=h_parms[1], dll_b=0.01,
                                       acquisition_T=acquisition_T, 
                                       fine_acq_T=fine_acq_T,
                                       tracking_T=tracking_T,
                                       M=M)
    ax.plot(trackresults.t, real.(trackresults.ZP))
    ax.plot(trackresults.t, imag.(trackresults.ZP))
end
subplots_adjust(hspace=0.4, wspace=0.4,  top=0.93, left=0.08, right=0.93, 
                bottom=0.01)
plot_title = string("Course Acq. T=$(floor(Int, acquisition_T*1000))ms;    ",
                    "Fine Acq. T=$(floor(Int, fine_acq_T*1000))ms;    ",
                    "Tracking T=$(floor(Int, tracking_T*1000))ms;    ",
                    "Noncoherent Integrations=$(M)")
suptitle(plot_title)


acqresults, trackresults, corr_result, SNR_est, 
P_d, above_threshold = process(data, 
                               sigtype, 
                               6, "Q";
                               fine_acq_method=:fft, 
                               return_corrresult=true,
                               return_Pd=true,
                               state_num=2, 
                               dynamickf=true,
                            #    h₀=h_parms[3], 
                            #    h₋₂=h_parms[1],
                               q_a=1,
                               dll_b=0.1,
                               acquisition_T=1e-3, 
                               fine_acq_T=100e-3,
                               tracking_T=20e-3,
                               M=20,
                               show_plot=true,
                               use_fine_acq=true);

#------------------------------------------------------------------------------#