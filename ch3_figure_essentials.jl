using GNSSTools
using Statistics
using PyPlot
pygui(true)


t = Array(0:0.1:100)
carrier = sin.((2π/10).*t)

# I channel
ranging_code_I = l5i_codes[1][1:10]
ranging_code_I[ranging_code_I .== 0] .= -1
neuman_code_I = nh10[5:9]
neuman_code_I[neuman_code_I .== 0] .= -1
navbits_I = [1, 1, -1]

# Q channel
ranging_code_Q = l5q_codes[1][1:10]
ranging_code_Q[ranging_code_Q .== 0] .= -1
neuman_code_Q = nh20[4:7]
neuman_code_Q[neuman_code_Q .== 0] .= -1
navbits_Q = [-1, -1, 1]

t_length = 1e-3
signal = definesignal(define_l1ca_code_type(), 5e6, 2*t_length; f_if=1.25e6,
                      include_thermal_noise=false, include_phase_noise=false,
                      include_adc=false, include_carrier=true)
generatesignal!(signal)
final_signal = real.(deepcopy(signal.data[1:200]))
final_signal = final_signal ./ maximum(final_signal)

signal = definesignal(define_l1ca_code_type(), 5e6, 2*t_length; f_if=1.25e6,
                      include_thermal_noise=true, include_phase_noise=true,
                      include_adc=false, include_carrier=true)
generatesignal!(signal)
final_signal_noise = real.(deepcopy(signal.data[1:200]))
final_signal_noise = final_signal_noise ./ maximum(final_signal_noise)

ranging_code_multiple = [ranging_code_I; ranging_code_I]

ylims = [-1.50, 1.50]
fig = figure(figsize=(6,3))
ax1 = fig.add_subplot(4,2,1)
ax1.step(ranging_code_I, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax2 = fig.add_subplot(4,2,2)
ax2.step(ranging_code_multiple, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax3 = fig.add_subplot(4,2,3)
ax3.plot(final_signal[1:100], "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax4 = fig.add_subplot(4,2,4)
ax4.plot(final_signal[101:200], "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax5 = fig.add_subplot(4,2,5)
ax5.plot(final_signal, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax6 = fig.add_subplot(4,2,6)
ax6.plot(final_signal_noise[1:100], "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax7 = fig.add_subplot(4,2,7)
ax7.plot(final_signal_noise[101:200], "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax8 = fig.add_subplot(4,2,8)
ax8.plot(final_signal_noise, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
# ax4 = fig.add_subplot(4,2,2)
# ax4.step(ranging_code_Q, "k")
# ylim(ylims)
# tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
# ax5 = fig.add_subplot(4,2,4)
# ax5.step(neuman_code_Q, "k")
# ylim(ylims)
# tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
# subplots_adjust(hspace=0.3, wspace=0.4, top=0.91, left=0.08, right=0.93, bottom=0.1)
savefig("figures/ch3_figure_essentials_signals.svg", dpi=300)


n0 = 3125
n0_offset = -25
n0_range = 100
signal = definesignal(define_l1ca_code_type(), 5e6, 1; f_if=1.25e6,
                      include_thermal_noise=true, include_phase_noise=true,
                      include_adc=false, include_carrier=true, 
                      code_start_idx=n0, CN0=40)
generatesignal!(signal)
replica = definesignal(define_l1ca_code_type(), 5e6, 1e-3)
fd_est, n0_est, SNR_est, corr_result1 = courseacquisition(signal, replica, 1; 
                                                         return_corrresult=true,
                                                         M=1, 
                                                         start_idx=1)
fd_est, n0_est, SNR_est, corr_result2 = courseacquisition(signal, replica, 1; 
                                                         return_corrresult=true,
                                                         M=1, 
                                                         start_idx=1+5000)
fd_est, n0_est, SNR_est, corr_result3 = courseacquisition(signal, replica, 1; 
                                                         return_corrresult=true,
                                                         M=10, 
                                                         start_idx=1)

corr_crop1 = corr_result1[:,n0-n0_range+n0_offset:n0+n0_range+n0_offset]
corr_crop2 = corr_result2[:,n0-n0_range+n0_offset:n0+n0_range+n0_offset]
corr_crop3 = corr_result3[:,n0-n0_range+n0_offset:n0+n0_range+n0_offset]

mplot3d = PyPlot.PyObject(PyPlot.axes3D)
fig = figure(figsize=(4,5))
ax1 = fig.add_subplot(3, 1, 1, projection="3d")
x, y = meshgrid(Array(1:size(corr_crop1)[2]), Array(1:size(corr_crop1)[1]))
surf(x, y, corr_crop1, antialiased=true, shade=true, alpha=0.9, cmap="coolwarm")
# ax1.grid(false)
# ax1.set_xticks([])
# ax1.set_yticks([])
# ax1.set_zticks([])
axis("off")
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)

ax3 = fig.add_subplot(3, 1, 2, projection="3d")
x, y = meshgrid(Array(1:size(corr_crop2)[2]), Array(1:size(corr_crop2)[1]))
surf(x, y, corr_crop2, antialiased=true, shade=true, alpha=0.9, cmap="coolwarm")
axis("off")
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)

ax3 = fig.add_subplot(3, 1, 3, projection="3d")
x, y = meshgrid(Array(1:size(corr_crop3)[2]), Array(1:size(corr_crop3)[1]))
surf(x, y, corr_crop3, antialiased=true, shade=true, alpha=0.9, cmap="coolwarm")
axis("off")
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)

savefig("figures/ch3_figure_essentials_acquisition_noncoherent.svg", dpi=300)


fd_est, n0_est, SNR_est, corr_result4 = courseacquisition(signal, replica, 1; 
                                                         return_corrresult=true,
                                                         M=2, 
                                                         start_idx=1)
replica2 = definesignal(define_l1ca_code_type(), 5e6, 2e-3)
fd_est, n0_est, SNR_est, corr_result5 = courseacquisition(signal, replica2, 1; 
                                                         return_corrresult=true,
                                                         M=1, 
                                                         start_idx=1)

fd_est, n0_est, SNR_est, corr_result6 = courseacquisition(signal, replica, 1; 
                                                         return_corrresult=true,
                                                         M=10, 
                                                         start_idx=1)
replica2 = definesignal(define_l1ca_code_type(), 5e6, 10e-3)
fd_est, n0_est, SNR_est, corr_result7 = courseacquisition(signal, replica2, 1; 
                                                         return_corrresult=true,
                                                         M=1, 
                                                         start_idx=1)

ylims = [-0.50, 1.50]
fig = figure(figsize=(6,3))
ax1 = fig.add_subplot(4,2,1)
result = corr_result1[argmax(corr_result1)[1],:]
ax1.plot(result ./ maximum(result), "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax2 = fig.add_subplot(4,2,2)
result = corr_result2[argmax(corr_result2)[1],:]
ax2.plot(result ./ maximum(result), "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax3 = fig.add_subplot(4,2,3)
result = corr_result5[argmax(corr_result5)[1],1:5000]
ax3.plot(result ./ maximum(result), "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax4 = fig.add_subplot(4,2,4)
result = corr_result4[argmax(corr_result4)[1],:]
ax4.plot(result ./ maximum(result), "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax5 = fig.add_subplot(4,2,5)
result = corr_result7[argmax(corr_result7)[1],1:5000]
ax5.plot(result ./ maximum(result), "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax6 = fig.add_subplot(4,2,6)
result = corr_result6[argmax(corr_result6)[1],:]
ax6.plot(result ./ maximum(result), "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
savefig("figures/ch3_figure_essentials_signals_correlation_results.svg", dpi=300)