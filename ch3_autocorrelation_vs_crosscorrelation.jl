using GNSSTools
using Statistics
using PyPlot
pygui(true)

N = 10
# idx = (rand(1:N), rand(1:N))
idx = (4, 4)
background_val = 1 
peak_val = 5
corr = fill(background_val, N, N)
corr[idx[1],idx[2]] = peak_val 

# Get codes for plot
unknown_code = circshift(deepcopy(l1ca_codes[1]), 500)
prn_1_code = deepcopy(l1ca_codes[1])
prn_2_code = deepcopy(l1ca_codes[2])
# Make zeros -1
A = 2
unknown_code[unknown_code .== 0] .= -A
prn_1_code[prn_1_code .== 0] .= -A
prn_2_code[prn_2_code .== 0] .= -A
unknown_code[unknown_code .== 1] .= A
prn_1_code[prn_1_code .== 1] .= A
prn_2_code[prn_2_code .== 1] .= A
# Perform auto and cross correlations
auto_correlation = real.(fft_correlate(unknown_code, prn_1_code))
cross_correlation = real.(fft_correlate(unknown_code, prn_2_code))

xlims = [-50, length(auto_correlation)+50]
ylims = [mean(cross_correlation)-120, maximum(auto_correlation)+53]

# # Simulate L1 C/A signal and perform course acquisition
# signal = definesignal(define_l1ca_code_type(), 5e6, 1; prn=1, f_d=1000,
#                       code_start_idx=1000, include_phase_noise=true)
# generatesignal!(signal)
# replica = definesignal(define_l1ca_code_type(), 5e6, 1e-3)
# fd_est, n0_est, SNR_est, corr_result = courseacquisition(signal, replica, 1; 
#                                                          return_corrresult=true,
#                                                          Δfd=1000)
# corr_crop = corr_result[:,n0_est-5:n0_est+5]
# idx_crop = argmax(corr_crop)
# box_x = [idx_crop[2]-0.5, idx_crop[2]-0.5, idx_crop[2]+0.5, idx_crop[2]+0.5, idx_crop[2]-0.5] .- 1
# box_y = [idx_crop[1]+0.5, idx_crop[1]-0.5, idx_crop[1]-0.5, idx_crop[1]+0.5, idx_crop[1]+0.5] .- 1

title_font_size = 9
axis_label_font_size = 8
ytick_vals = [0, 1023]
xtick_vals = [1, 501, 1023]
xtick_labels = [string(i) for i in [524, 1, 523]]
                                                         
fig = figure(figsize=(6.5,1.5))
ax1 = fig.add_subplot(1,2,1)
ax1.plot(auto_correlation, "k")
xlim(xlims)
ylim(ylims)
ax1.set_yticks(ytick_vals)
ax1.set_xticks(xtick_vals)
ax1.set_xticklabels(xtick_labels)
yticks(fontsize=axis_label_font_size)
xticks(fontsize=axis_label_font_size)
xlabel("Chips", fontsize=axis_label_font_size)
ylabel("rⁱ[n]  ⋆  sⁱ[n]", fontsize=axis_label_font_size)
title("(a)", fontsize=title_font_size)
# xlabel("Time", fontsize=axis_label_font_size)
tick_params(left=true, bottom=true, labelleft=true, labelbottom=true)
ax2 = fig.add_subplot(1,2,2)
ax2.plot(cross_correlation, "k")
xlim(xlims)
ylim(ylims)
ax2.set_yticks(ytick_vals)
ax2.set_xticks(xtick_vals)
ax2.set_xticklabels(xtick_labels)
yticks(fontsize=axis_label_font_size)
xticks(fontsize=axis_label_font_size)
xlabel("Chips", fontsize=axis_label_font_size)
ylabel("rⁱ[n]  ⋆  sʲ[n]", fontsize=axis_label_font_size)
title("(b)", fontsize=title_font_size)
tick_params(left=true, bottom=true, labelleft=true, labelbottom=true)
# tick_params(left=true, bottom=false, labelleft=true, labelbottom=false)
# tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
subplots_adjust(hspace=0.5, wspace=0.4, top=0.85, left=0.1, right=0.93, bottom=0.3)
savefig("figures/ch3_autocorrelation_vs_cross_correlation.svg", dpi=300)


# Simulate L1 C/A signal and perform course acquisition
signal = definesignal(define_l1ca_code_type(), 5e6, 1; prn=1, f_d=1000,
                      code_start_idx=1000, include_phase_noise=true)
generatesignal!(signal)
replica = definesignal(define_l1ca_code_type(), 5e6, 1; include_phase_noise=false,
                       include_thermal_noise=false, include_adc=false,
                       )