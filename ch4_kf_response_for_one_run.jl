using GNSSTools
using ProgressMeter
using JLD
using PyPlot
pygui(true)
include("kalman_filter_gain_error_functions.jl")


function kf_run_results(sigtype, f_s, t_length;
                        f_d=0, fd_rate=0, phi_init=0, f_if=0,
                        include_phase_noise=true, include_thermal_noise=true,
                        include_adc=true, include_carrier=true,
                        acquisition_T=1e-3, M=1, fine_acq_T=10e-3,
                        tracking_T=1e-3, fine_acq_method=:fft,
                        CN0=45, code_start_idx=1, fd_range=5000, prn=1,
                        σᵩ²=GNSSTools.phase_noise_variance(CN0, tracking_T),
                        h_parms=[1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2,
                        fd_center=0, nADC=8, dynamickf=true, state_num=3,
                        Tsys=535, channel="I", q_a=1, use_fine_acq=true,
                        dll_b=0.1, p_d=0.9, p_fa=1e-7, σω=1,
                        show_process_result_plot=false, fd_rate_init=0,
                        figsize=(8, 7.5), window_size=10, 
                        save_to_file=missing, top=0.9, bottom=0.1,
                        left=0.1, right=0.9, wspace=0.15, hspace=0.4)
    sig_freq = sigtype.sig_freq
    signal = definesignal(sigtype, f_s, t_length; f_d=f_d,
                          prn=prn, fd_rate=fd_rate, f_if=f_if,
                          include_phase_noise=include_phase_noise, 
                          include_thermal_noise=include_thermal_noise, 
                          include_adc=include_adc, receiver_h_parms=h_parms,
                          include_carrier=include_carrier, nADC=nADC,
                          CN0=CN0, phi=phi_init, Tsys=Tsys, 
                          code_start_idx=code_start_idx)
    generatesignal!(signal)
    results, trackresults, corr_result, SNR_est, 
    P_d, above_threshold = process(signal, sigtype, 1, channel;
                                   return_corrresult=true,
                                   return_Pd=true,
                                   show_plot=show_process_result_plot,
                                   acquisition_T=acquisition_T,
                                   M=M,
                                   fine_acq_method=fine_acq_method,
                                   use_fine_acq=use_fine_acq,
                                   fine_acq_T=fine_acq_T,
                                   tracking_T=tracking_T,
                                   h₋₂=h₋₂,
                                   h₀=h₀,
                                   σᵩ²=σᵩ²,
                                   q_a=q_a,
                                   state_num=state_num,
                                   dynamickf=dynamickf,
                                   fd_rate=fd_rate_init,
                                   dll_b=dll_b,
                                   p_d=p_d,
                                   p_fa=p_fa,
                                   σω=σω)
    
    # Get the matrices used in signal tracking 
    A = trackresults.A
    C = trackresults.C
    Q = trackresults.Q
    R = trackresults.R
    # Get KF gain values
    t = trackresults.t
    ZP = trackresults.ZP
    K = trackresults.K
    k₁ = K[1,:]
    k₂ = K[2,:]
    if state_num == 3
        k₃ = K[3,:]
    end
    B_exp_numerical = zeros(length(t))
    B_exp_analytical = zeros(length(t))
    t_window = t[window_size:window_size:end]
    dphi_meas = trackresults.dphi_meas
    dphi_filt = trackresults.dphi_filt
    σᵩ_unfilt = zeros(length(t_window))
    σᵩ_filt = zeros(length(t_window))
    pᵩ = zeros(length(t))
    for i in 1:length(t)
        # Calculate experimental steady state gain using the numerical solution
        B_exp_numerical[i] = get_B(init_H(K[:,i]...))
        # Calculate experimental steady state gain using the analytical equation
        B_exp_analytical[i] = calc_B(K[:,i]...)
        # Calculate experimental pᵩ
        P = trackresults.P_full[i,:,:]
        pᵩ[i] = sqrt((C*P*C')[1])
        if i%window_size == 0
            # Calculate standard deviation of unfiltered carrier phase error
            σᵩ_unfilt[Int(i/window_size)] = std(dphi_meas[i-(window_size-1):i])
            # Calculate standard deviation of filtered carrier phase error
            σᵩ_filt[Int(i/window_size)] = std(dphi_filt[i-(window_size-1):i])
        end
    end

    # Calculate the theoretical steady state kalman gain
    K_steady_state = dkalman(A, C, Q, GNSSTools.Diagonal(R))
    k₁_steady_state = K_steady_state[1,:]
    k₂_steady_state = K_steady_state[2,:]
    if state_num == 3
        k₃_steady_state = K_steady_state[3,:]
    end
    # Calculate theoretical steady state gain using the numerical solution
    B_theo_numerical = get_B(init_H(K_steady_state...))
    # Calculate theoretical steady state gain using the analytical equation
    B_theo_analytical = calc_B(K_steady_state...)
    # Calculate theoretical pᵩ
    P_steady_state = GNSSTools.dare(A', C', Q, GNSSTools.Diagonal(R)')
    pᵩ_steady_state = sqrt((C*P_steady_state*C')[1])
    # Calculate the steady state unfiltered carrier phase error (σᵩ)
    σᵩ_unfilt_steady_state = sqrt(R[1])
    σᵩ_filt_steady_state_num = NaN
    σᵩ_filt_steady_state_ana = NaN
    if include_phase_noise  # Calculate σᵩ with oscillator phase noise
        # Calculate the steady state filtered carrier phase error using
        # numerical solution
        σᵩ_filt_steady_state_num = sqrt(phase_noise_variance(B_theo_numerical, 
                                                            CN0, tracking_T, 
                                                            h₀, h₋₁, h₋₂))
        # Calculate the steady state filtered carrier phase error using 
        # analytical equation
        try
            σᵩ_filt_steady_state_ana = sqrt(phase_noise_variance(B_theo_analytical, 
                                                                CN0, tracking_T, 
                                                                h₀, h₋₁, h₋₂))
        catch
            σᵩ_filt_steady_state_ana = NaN
            @warn("Unable to calculate value for `σᵩ_filt_steady_state_ana` due to negative sqrt.")
        end
    else # Calculate σᵩ without oscillator phase noise
        # Calculate the steady state filtered carrier phase error using
        # numerical solution
        σᵩ_filt_steady_state_num = sqrt(phase_noise_variance(B_theo_numerical, 
                                                            CN0, tracking_T))
        # Calculate the steady state filtered carrier phase error using 
        # analytical equation
        try
            σᵩ_filt_steady_state_ana = sqrt(phase_noise_variance(B_theo_analytical, 
                                                                CN0, tracking_T))
        catch
            σᵩ_filt_steady_state_ana = NaN
            @warn("Unable to calculate value for `σᵩ_filt_steady_state_ana` due to negative sqrt.")
        end
    end

    # Make plot
    fig = figure(figsize=figsize)
    rows = 3
    cols = 2
    
    # KF gains
    ax1 = fig.add_subplot(rows, cols, 1)
    ax1.plot(t, k₁, "-", color="blue", label="k₁")
    ax1.plot(t, k₂, "-", color="orange", label="k₂")
    ax1.plot([t[1], t[end]], fill(k₁_steady_state, 2), ":", color="blue")
    ax1.plot([t[1], t[end]], fill(k₂_steady_state, 2), ":", color="orange")
    if state_num == 3
        ax1.plot(t, k₃, "-", color="green", label="k₃")
        ax1.plot([t[1], t[end]], fill(k₃_steady_state, 2), ":", color="green")
    end
    ax1.plot([], [], "k:", label="Steady State kₙ")
    xlabel("Time (s)")
    ylabel("KF Gain")
    legend()
    title("(a)")

    # Filtered carrier phase errors based off P
    ax2 = fig.add_subplot(rows, cols, 2)
    ax2.plot(t, pᵩ.*(180/π), "k-", label="pᵩ")
    ax2.plot([t[1], t[end]], fill(pᵩ_steady_state, 2).*(180/π), "k:", 
             label="Steady State pᵩ")
    xlabel("Time (s)")
    ylabel("pᵩ (degrees)")
    legend()
    title("(b)")

    # Unfiltered carrier phase error
    ax3 = fig.add_subplot(rows, cols, 3)
    ax3.plot(t_window, σᵩ_unfilt.*(180/π), "k-", label="σᵩ Meas.")
    ax3.plot([t[1], t[end]], fill(σᵩ_unfilt_steady_state, 2).*(180/π), "k:", 
              label="σᵩ Theo.")
    xlabel("Time (s)")
    ylabel("Unfiltered σᵩ (degrees)")
    legend()
    title("(c)")

    # Filtered carrier phase error
    ax4 = fig.add_subplot(rows, cols, 4)
    ax4.plot(t_window, σᵩ_filt.*(180/π), "k-", label="σᵩ Meas.")
    ax4.plot([t[1], t[end]], fill(σᵩ_filt_steady_state_num, 2).*(180/π), "k:", 
              label="σᵩ Theo. Num.")
    # ax4.plot([t[1], t[end]], fill(σᵩ_filt_steady_state_ana, 2).*(180/π), ":", 
            #   label="σᵩ Theo. Ana.")
    xlabel("Time (s)")
    ylabel("Filtered σᵩ (degrees)")
    legend()
    title("(d)")

    # Estimated equivalent KF noise bandwidth
    ax5 = fig.add_subplot(rows, cols, 5)
    ax5.plot(t, B_exp_numerical, "k-", label="B Exp. Num.")
    # ax5.plot(t, B_exp_analytical, ":", color="blue", label="B Exp. Ana.")
    ax5.plot([t[1], t[end]], fill(B_theo_numerical, 2), "k:", 
             label="B Theo. Num.")
    # ax5.plot([t[1], t[end]], fill(B_theo_analytical, 2), ":", color="orange",
            #  label="B Theo. Ana.")
    xlabel("Time (s)")
    ylabel("Bandwidth (Hz)")
    ylim([-2, 15])
    legend()
    title("(e)")

    # Prompt correlator outputs
    ax6 = fig.add_subplot(rows, cols, 6)
    ax6.plot(t, real.(ZP), label="real(Zₚ)")
    ax6.plot(t, imag.(ZP), label="imag(Zₚ)")
    xlabel("Time (s)")
    ylabel("Zₚ")
    legend()
    title("(f)")

    subplots_adjust(hspace=hspace, wspace=wspace,  top=top, left=left, right=right, 
                    bottom=bottom)
    if ~ismissing(save_to_file)
        savefig(save_to_file, dpi=300)
    end
end


tcxo_h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15]./16.3676e6^2  # TCXO
ocxo_h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
h_parms = tcxo_h_parms

t_length = 1
f_s = 2.5e6
f_if = 0
f_d = 1000
fd_rate = 0
phi_init = 0
code_start_idx = 1
CN0 = 45
nADC = 8

sigtype = define_l1ca_code_type()
channel = "I"

dll_b = 1
acquisition_T = 40e-3
M = 5
fine_acq_T = 100e-3
fine_acq_method = :fft
tracking_T = 1e-3
q_a = 0.1
σω = 1
state_num = 2
window_size = 50
fd_rate_init = 0

include_thermal_noise = true
include_phase_noise = true
include_adc = true
include_carrier = true

show_process_result_plot = false

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/figures/"
save_to_file = string(directory, "ch4_kf_response_successful_track_$(state_num)-state.pdf")
# save_to_file = string(directory, "ch4_kf_response_unsuccessful_track_$(state_num)-state.pdf")

kf_run_results(sigtype, f_s, t_length;
                        f_d=f_d, fd_rate=fd_rate, phi_init=phi_init, f_if=f_if,
                        include_phase_noise=include_phase_noise, 
                        include_thermal_noise=include_thermal_noise,
                        include_adc=include_adc, include_carrier=include_carrier,
                        acquisition_T=acquisition_T, M=M, 
                        fine_acq_T=fine_acq_T,
                        tracking_T=tracking_T, 
                        fine_acq_method=fine_acq_method,
                        CN0=CN0, code_start_idx=code_start_idx,
                        h_parms=h_parms, nADC=nADC, 
                        state_num=state_num,
                        channel=channel, q_a=q_a,
                        dll_b=dll_b,
                        show_process_result_plot=show_process_result_plot,
                        figsize=(8, 7.5), window_size=window_size, 
                        top=0.95, bottom=0.1, left=0.1, right=0.9, 
                        wspace=0.25, hspace=0.52, fd_rate_init=fd_rate_init,
                        save_to_file=save_to_file)
