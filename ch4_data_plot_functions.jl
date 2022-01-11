function plot_example_track_results(trackresults; f_if=0, M=10,
                                    save_to_file=missing,
                                    code_length=1023,
                                    top=0.9,
                                    bottom=0.1,
                                    left=0.1,
                                    right=0.9,
                                    wspace=0.1,
                                    hspace=0.1,
                                    marker_size=1,
                                    line_width=1,
                                    use_ylims=true,
                                    truth_snr=missing,
                                    truth_doppler=missing,
                                    truth_doppler_rate=missing,
                                    truth_t=missing,
                                    truth_phi=missing,
                                    truth_code_phase=missing)
    t = trackresults.t
    ZP = trackresults.ZP
    SNR = trackresults.SNR
    x = trackresults.x
    P = trackresults.P
    ϕ = x[1,:] ./ 2π #.* (180/π)
    σ_ϕ = sqrt.(P[1,:]) ./ 2π # .* (180/π)
    f_d = x[2,:] ./ (2π) .- f_if
    σ_f_d = sqrt.(P[2,:]) ./ (2π)
    f_dd = x[3,:] ./ (2π)
    σ_f_dd = sqrt.(P[3,:]) ./ (2π)
    dphi_meas = trackresults.dphi_meas
    dphi_filt = trackresults.dphi_filt
    code_err_meas = trackresults.code_err_meas
    code_err_filt = trackresults.code_err_filt
    code_phase_meas = trackresults.code_phase_meas
    code_phase_filt = trackresults.code_phase_filt
    N_over_2 = floor(Int, length(t)/2)
    fig = figure(figsize=(8, 7.5))
    ax1 = fig.add_subplot(4,2,1)
    ax1.plot(t, real.(ZP), label="real(ZP)", linewidth=line_width)
    ax1.plot(t, imag.(ZP), label="imag(ZP)", linewidth=line_width)
    xlabel("Time (seconds)")
    ylabel("ZP")
    legend()
    title("(a)")
    ax2 = fig.add_subplot(4,2,2)
    if ~ismissing(truth_snr)
        ax2.plot(t, SNR, "k.", markersize=marker_size, label="Measured SNR")
        ax2.axhline(y=truth_snr, c="blue", linestyle="-", label="Truth SNR",
                    linewidth=line_width)
        legend()
    else
        ax2.plot(t, SNR, "k.", markersize=marker_size)
    end
    xlabel("Time (seconds)")
    ylabel("SNR (dB)")
    title("(b)")
    ax3 = fig.add_subplot(4,2,3)
    if ~ismissing(truth_phi)
        phi_diff = truth_phi./2π .- ϕ
        phi_diff = phi_diff .- mean(phi_diff)
        ax3.plot(t,  phi_diff, "k-", label=L"ϕ - \hat{\phi}", linewidth=line_width)
        ax3.plot(t, -σ_ϕ, label=L"\pm\sigma_{\phi}",
                 color="grey", linestyle=":", linewidth=line_width)
        ax3.plot(t, σ_ϕ, color="grey", linestyle=":", linewidth=line_width)
        ylabel(string(L"ϕ - \hat{\phi}", " (cycles)"))
    else
        ax3.plot(t, ϕ, "k-", label=L"\hat{\phi}", linewidth=line_width)
        ax3.plot(t, (ϕ .- σ_ϕ), label=L"\pm\sigma_{\phi}",
                 color="grey", linestyle=":", linewidth=line_width)
        ax3.plot(t, (ϕ .+ σ_ϕ), color="grey", linestyle=":", linewidth=line_width)
        ylabel(string(L"\hat{\phi}", " (cycles)"))
    end
    xlabel("Time (seconds)")
    legend()
    title("(c)")
    ax4 = fig.add_subplot(4,2,4)
    ax4.plot(t, dphi_meas .* (180/π), "k.", label=L"\Delta\phi_{measured}", markersize=marker_size)
    ax4.plot(t, dphi_filt .* (180/π), "b-", label=L"\Delta\phi_{filtered}", linewidth=line_width)
    xlabel("Time (seconds)")
    ylabel(string(L"\Delta\phi", " (degrees)"))
    legend()
    title("(d)")
    ax5 = fig.add_subplot(4,2,5)
    ax5.plot(t, f_d, "k-", label=L"\hat{f}_d", linewidth=line_width)
    ax5.plot(t, f_d .- σ_f_d, label=L"\pm\sigma_{f_d}",
             color="grey", linestyle=":", linewidth=line_width)
    ax5.plot(t, f_d .+ σ_f_d, color="grey", linestyle=":", linewidth=line_width)
    if ~ismissing(truth_doppler)
        ax5.plot(truth_t, truth_doppler, "b-", linewidth=line_width, 
                 label=string("Truth ", L"f_d"))
    end
    if use_ylims
        ylim([median(f_d .- σ_f_d.*8), median(f_d .+ σ_f_d.*8)])
    end
    xlabel("Time (seconds)")
    ylabel(string(L"\hat{f}_d", " (Hz)"))
    legend()
    title("(e)")
    ax6 = fig.add_subplot(4,2,6)
    ax6.plot(t, code_err_meas, "k.", label=L"\Delta n_{measured}", markersize=marker_size)
    ax6.plot(t, code_err_filt, "b-", label=L"\Delta n_{filtered}", linewidth=line_width)
    xlabel("Time (seconds)")
    ylabel("Code Phase Error (chips)")
    legend()
    title("(f)")
    ax7 = fig.add_subplot(4,2,7)
    ax7.plot(t, f_dd, "k-", label=L"\hat{\dot{f}}_d", linewidth=line_width)
    ax7.plot(t, f_dd .- σ_f_dd, label=L"\pm\sigma_{\dot{f}_d}",
             color="grey", linestyle=":", linewidth=line_width)
    ax7.plot(t, f_dd .+ σ_f_dd, color="grey", linestyle=":", linewidth=line_width)
    if ~ismissing(truth_doppler_rate)
        ax7.plot(truth_t, truth_doppler_rate, "b-", linewidth=line_width, 
                 label=string("Truth ", L"\dot{f}_d"))
    end
    if use_ylims
        ylim([median(f_dd .- σ_f_dd.*4), median(f_dd .+ σ_f_dd.*4)])
    end
    # ylim([-408, 88])
    xlabel("Time (seconds)")
    ylabel(string(L"\hat{\dot{f}}_d", " (Hz/s)"))
    legend()
    title("(g)")
    ax8 = fig.add_subplot(4,2,8)
    ax8.plot(t, code_phase_filt .% code_length, "k-", label=L"\hat{n}", linewidth=line_width)
    intervals = floor(Int, length(t)/M)
    t_sigma_n = Array{Float64}(undef, intervals+1)
    sigma_n = Array{Float64}(undef, intervals+1)
    n = Array{Float64}(undef, intervals+1)
    for i in 1:intervals
        sigma_n[i+1] = std(code_err_meas[(i-1)*M+1:i*M])
        # n[i] = mean(code_phase_filt[(i-1)*M+1:i*M])
        # n[i] = code_phase_filt[(i-1)*M+floor(Int, M/2)]
        n[i+1] = code_phase_filt[i*M]
        t_sigma_n[i+1] = t[i*M]
        if i == 1
            t_sigma_n[1] = t[1]
            n[1] = code_phase_filt[1]
            sigma_n[1] = sigma_n[i+1]
        end
    end
    ax8.plot(t_sigma_n, (n .- sigma_n) .% code_length, label=L"\pm\sigma_n",
             color="grey", linestyle=":", linewidth=line_width)
    ax8.plot(t_sigma_n, (n .+ sigma_n) .% code_length, color="grey", linestyle=":", 
             linewidth=line_width)
    if ~ismissing(truth_code_phase)
        ax8.plot(truth_t, truth_code_phase .% code_length, "b-", linewidth=line_width, 
                 label=string("Truth ", L"n"))
    end
    xlabel("Time (seconds)")
    ylabel("Code Phase (chips)")
    legend()
    title("(h)")
    subplots_adjust(hspace=hspace, wspace=wspace,  top=top, left=left, right=right, 
                    bottom=bottom)
    if ~ismissing(save_to_file)
        savefig(save_to_file, dpi=300)
    end
end