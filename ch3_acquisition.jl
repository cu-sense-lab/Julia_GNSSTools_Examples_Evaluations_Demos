using GNSSTools
using Statistics
using PyPlot
pygui(true)


CN0 = 36
prn = 10
t_length = 1
fd_range = 5000
f_s = 5e6
sigtype = define_l1ca_code_type()
# Simulate L1 C/A signal and perform course acquisition
signal = definesignal(sigtype, f_s, t_length; prn=prn, f_d=4910,
                      code_start_idx=2500, include_phase_noise=true, CN0=CN0)
generatesignal!(signal)
replica = definesignal(sigtype, f_s, 1e-3)
Δfd = 1/replica.t_length
fd_est, n0_est, SNR_est, corr_result = courseacquisition(signal, replica, prn; 
                                                         return_corrresult=true,
                                                         Δfd=Δfd, fd_range=fd_range,
                                                         M=1)

max_idx = argmax(corr_result)
ps = corr_result[max_idx[1],:]
PS = corr_result[max_idx]
N, M = size(corr_result)
noise_val = sum(corr_result[max_idx[1],:])
PN = (noise_val - (replica.t_length/1e-3)*PS)/(M - 1)
SNR = 10log10(PS/PN)

################################################################################
h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
file = "/home/sjbilardi/Documents/SeNSeLab/GNSSTools.jl/test/data/known_oscillator_data/20210327_013800_000005_Aero-North-TCXO_1575.42_25.sc16"
data = loaddata(:sc16, file, 25e6, L1_freq, L1_freq, 1., skip_to=1+0.00031216)


replica = definesignal(define_l1ca_code_type(), 25e6, 2e-3)
for i in 1:32
    fd_est, n0_est, SNR_est = courseacquisition(data, replica, i; M=10)
    println("PRN $i\t$SNR_est")
end

fig = figure()
acquisiton_T = 1e-3
fine_acq_T = 10e-3
tracking_T = 1e-3
for i in 1:32
    ax = fig.add_subplot(8, 4, i)
    tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    title("PRN $i")
    acqresults, trackresults = process(data, define_l1ca_code_type(), i;
                                       show_plot=false,
                                       fine_acq_method=:fft, state_num=3,
                                       dynamickf=true, h₀=h_parms[3], 
                                       h₋₂=h_parms[1], dll_b=1,
                                       acquisition_T=acquisiton_T, 
                                       fine_acq_T=fine_acq_T,
                                       tracking_T=tracking_T,
                                       M=10)
    ax.plot(trackresults.t, real.(trackresults.ZP))
    ax.plot(trackresults.t, imag.(trackresults.ZP))
end
subplots_adjust(hspace=0.4, wspace=0.4,  top=0.93, left=0.08, right=0.93, 
                bottom=0.01)
plot_title = string("Course Acq. T=$(floor(Int, acquisiton_T*1000))ms;    ",
                    "Fine Acq. T=$(floor(Int, fine_acq_T*1000))ms;    ",
                    "Tracking T=$(floor(Int, tracking_T*1000))ms")
suptitle(plot_title)