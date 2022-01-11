using GNSSTools
using Statistics
using PyPlot
pygui(true)
include("ch4_data_plot_functions.jl")

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
folder = "figures/l1ca_data/"

file = "/home/sjbilardi/Documents/SeNSeLab/GNSSTools.jl/test/data/20201211_211745_GPSL1_1575.42_25.0.sc4"
data = loaddata(:sc4, file, 25e6, L1_freq, L1_freq, 2, skip_to=1)

replica = definereplica(define_l1ca_code_type(), 25e6, 2e-3)
detectable_prns = []
for i in 1:32
    fd_est, n0_est, SNR_est, 
    P_d, above_threshold = courseacquisition(data, replica, i; M=10, 
                                             Δfd=1/(1*replica.t_length))
    if above_threshold && (P_d >= 0.9)
        detectable = true
        push!(detectable_prns, i)
    else
        detectable = false
    end
    f_code_d, f_code_dd = GNSSTools.calc_doppler_code_rate(l1ca_chipping_rate, 
                                                           replica.signal_type.sig_freq, 
                                                           fd_est, 0)
    n0 = calcinitcodephase(l1ca_code_length, f_code_d, f_code_dd, data.f_s, 
                           n0_est)
    println("PRN $i\t$SNR_est\t$fd_est\t$n0_est\t$n0\t$(Float64(round(P_d, digits=2)))\t$detectable")
end
println("$(length(detectable_prns)) Detected PRNs: $detectable_prns")

fig = figure(figsize=(7.5, 5))
acquisition_T = 2e-3
fine_acq_T = 10e-3
fine_acq_sub_T = 1e-3
tracking_T = 1e-3
M = 10
linewidth = 0.5
tracking_σᵩ = Array{Float64}(undef, length(detectable_prns))
tracking_snr = Array{Float64}(undef, length(detectable_prns))
tracking_code_phase_err = Array{Float64}(undef, length(detectable_prns))
for i in 1:length(detectable_prns)
    prn = detectable_prns[i]
    ax = fig.add_subplot(3, 3, i)
    tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    title("PRN $prn")
    acqresults, trackresults = process(data, define_l1ca_code_type(), prn;
                                       show_plot=false,
                                       fine_acq_method=:fft, 
                                       state_num=2,
                                       dynamickf=true, 
                                    #    h₀=h_parms[3], h₋₂=h_parms[1], 
                                       dll_b=0.01,
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
savefig(string(directory, folder, "ch4_l1ca_acquired_prns.svg"), dpi=300)


for prn in detectable_prns
    acqresults, trackresults, corr_result, SNR_est, 
    P_d, above_threshold = process(data, 
                                define_l1ca_code_type(), 
                                prn;
                                fine_acq_method=:fft, 
                                return_corrresult=true,
                                return_Pd=true,
                                state_num=3, 
                                #    h₀=h_parms[3], 
                                #    h₋₂=h_parms[1],
                                q_a=1,
                                dll_b=1,
                                acquisition_T=2e-3, 
                                fine_acq_T=10e-3,
                                fine_acq_sub_T=1e-3,
                                tracking_T=1e-3,
                                show_plot=false,
                                σω=1000,
                                M=10,
                                use_fine_acq=true);

    name = "pdfs/ch4_l1ca_prn_$(trackresults.prn).pdf"
    file_name = string(directory, folder, name)
    plot_example_track_results(trackresults; f_if=0, code_length=l1ca_code_length,
                            save_to_file=file_name,
                            bottom=0.075, top=0.95, left=0.125, right=0.875, 
                            hspace=0.75, wspace=0.35, marker_size=1,
                            line_width=0.8)
end

