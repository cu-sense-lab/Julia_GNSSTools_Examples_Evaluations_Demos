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
signal = definesignal(define_l1ca_code_type(), 5e6, t_length; f_if=1.25e6,
                      include_thermal_noise=false, include_phase_noise=false,
                      include_adc=false, include_carrier=true)
generatesignal!(signal)
final_signal = real.(signal.data[1:100])
final_signal = final_signal ./ maximum(final_signal)

ylims = [-1.50, 1.50]
fig = figure(figsize=(6,3))
ax1 = fig.add_subplot(4,2,1)
ax1.step(ranging_code_I, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax2 = fig.add_subplot(4,2,3)
ax2.step(neuman_code_I, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax3 = fig.add_subplot(4,2,5)
ax3.step(navbits_I, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax4 = fig.add_subplot(4,2,2)
ax4.step(ranging_code_Q, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax5 = fig.add_subplot(4,2,4)
ax5.step(neuman_code_Q, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax6 = fig.add_subplot(4,2,6)
ax6.step(navbits_Q, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax7 = fig.add_subplot(4,2,7)
ax7.plot(carrier, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax8 = fig.add_subplot(4,2,8)
ax8.plot(final_signal, "k")
ylim(ylims)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
# subplots_adjust(hspace=0.3, wspace=0.4, top=0.91, left=0.08, right=0.93, bottom=0.1)
savefig("figures/ch2_custom_code_def_waveforms.svg", dpi=300)


