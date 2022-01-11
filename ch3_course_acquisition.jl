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
unknown_code[unknown_code .== 0] .= -1
prn_1_code[prn_1_code .== 0] .= -1
prn_2_code[prn_2_code .== 0] .= -1
# Perform auto and cross correlations
auto_correlation = real.(fft_correlate(unknown_code, prn_1_code))
cross_correlation = real.(fft_correlate(unknown_code, prn_2_code))

xlims = [-50, length(auto_correlation)+50]
ylims = [mean(cross_correlation)-120, maximum(auto_correlation)+53]

# Simulate L1 C/A signal and perform course acquisition
signal = definesignal(define_l1ca_code_type(), 5e6, 1; prn=1, f_d=1000,
                      code_start_idx=1000, include_phase_noise=true)
generatesignal!(signal)
replica = definesignal(define_l1ca_code_type(), 5e6, 1e-3)
fd_est, n0_est, SNR_est, corr_result = courseacquisition(signal, replica, 1; 
                                                         return_corrresult=true,
                                                         Δfd=1000)
corr_crop = corr_result[:,n0_est-5:n0_est+5]
idx_crop = argmax(corr_crop)
box_x = [idx_crop[2]-0.5, idx_crop[2]-0.5, idx_crop[2]+0.5, idx_crop[2]+0.5, idx_crop[2]-0.5] .- 1
box_y = [idx_crop[1]+0.5, idx_crop[1]-0.5, idx_crop[1]-0.5, idx_crop[1]+0.5, idx_crop[1]+0.5] .- 1

title_font_size = 11
axis_label_font_size = 10
                                                         
fig = figure(figsize=(6.5,2))
ax1 = subplot2grid((2,4), (0,0), rowspan=2, colspan=2)
ax1.imshow(corr_crop, aspect="equal", cmap="gray", vmin=-9*mean(corr_crop), 
           vmax=1.5*maximum(corr_crop))
ax1.plot(box_x, box_y, "r")
ax1.set_xticks(Array(1:N).-0.5)
ax1.set_yticks(Array(1:N).-0.5)
rc("grid", linestyle=":", color="gray", linewidth=1)
grid(true)
xlabel("Code Phase", fontsize=axis_label_font_size)
ylabel("Doppler Frequency", fontsize=axis_label_font_size)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
title("(a)")
ax2 = subplot2grid((2,4), (0,2), rowspan=1, colspan=2)
ax2.plot(auto_correlation, "k")
xlim(xlims)
ylim(ylims)
title("(b)")
# xlabel("Time", fontsize=axis_label_font_size)
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
ax3 = subplot2grid((2,4), (1,2), rowspan=1, colspan=2)
ax3.plot(cross_correlation, "k")
xlim(xlims)
ylim(ylims)
xlabel("Time", fontsize=axis_label_font_size)
title("(c)")
tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
subplots_adjust(hspace=0.5, wspace=0.4, top=0.85, left=0.08, right=0.93, bottom=0.1)
savefig("figures/ch3_course_acquisition.pdf", dpi=300)