using GNSSTools
using ProgressMeter
using JLD
using Statistics
using PyPlot
pygui(true)


function get_metric(trackresults; start_t=0.5)
    T = trackresults.T
    start_idx = floor(Int, start_t/T)
    σ²ᵩ = var(trackresults.dphi_meas[start_idx:end])
    return σ²ᵩ
end


function compare_tracking_cases(f_d_cases, fd_rate_cases, h_parms, CN0s, Ts,
                                sigtype, f_s, t_length,
                                qₐ, state_num, σω, acquisition_T, M,
                                fine_acq_T, dll_b; 
                                doppler_bin_num=2, f_if=0, nADC=8,
                                iterations=1, p_fa=1e-7, p_d=0.9, prn=1,
                                attempts=3, save_to_file=missing, channel="I",
                                case_names=missing,
                                use_fd_rate_init=false) 
    results = zeros(length(Ts), length(f_d_cases), length(h_parms), length(CN0s))
    signal = definesignal(sigtype, f_s, t_length; f_if=f_if, nADC=nADC)
    N = length(Ts)*length(f_d_cases)*length(h_parms)*length(CN0s)*iterations
    sample_num_per_1ms = floor(Int, f_s*1e-3)
    p = Progress(N, 1, "Processing...")
    for i in 1:length(Ts)
        T = Ts[i]
        for j in 1:length(f_d_cases)
            f_d = f_d_cases[j]
            fd_center = round(f_d*T)/T
            fd_range = doppler_bin_num/T
            fd_rate = fd_rate_cases[j]
            for k in 1:length(h_parms)
                h_parm = h_parms[1]
                h₀ = h_parm[3]
                h₋₁ = h_parm[2]
                h₋₂ = h_parm[1]
                for g in 1:length(CN0s)
                    CN0 = CN0s[g]
                    metric_val = 0.
                    for iteration in 1:iterations
                        P_d = 0
                        # attempt = 0
                        # while (P_d < p_d) || (attempt != attempts)
                            ϕ₀ = rand(0:0.0001:2π)
                            code_start_idx = rand(1:1:sample_num_per_1ms)
                            definesignal!(signal; f_d=f_d, fd_rate=fd_rate, 
                                          phi=ϕ₀, code_start_idx=code_start_idx, 
                                          receiver_h_parms=h_parm, CN0=CN0,
                                          new_thermal_noise=true, 
                                          new_phase_noise=true)
                            generatesignal!(signal)
                            σᵩ² = phase_noise_variance(CN0, T)
                            if use_fd_rate_init
                                fd_rate_init = fd_rate + σω*randn()
                            else
                                fd_rate_init = 0
                            end
                            acqresults, trackresults, P_d, above_threshold = process(signal, 
                                            sigtype, 
                                            prn, channel;
                                            fine_acq_method=:fft, 
                                            fd_center=fd_center,
                                            fd_range=fd_range, 
                                            h₀=h₀, 
                                            h₋₂=h₋₂, 
                                            q_mult=1,
                                            cov_mult=1,
                                            R_mult=1,
                                            σᵩ²=σᵩ²,
                                            σω=σω,
                                            q_a=qₐ,
                                            dynamickf=true,
                                            state_num=state_num,
                                            dll_b=dll_b,
                                            acquisition_T=acquisition_T,
                                            fine_acq_T=fine_acq_T,
                                            tracking_T=T,
                                            M=M,
                                            fd_rate=fd_rate_init,
                                            show_plot=false,
                                            return_Pd=true)
                            # attempt += 1
                        # end
                        next!(p)
                        metric_val += get_metric(trackresults)
                    end
                    results[i,j,k,g] = sqrt(metric_val / iterations)
                end
            end
        end
    end
    if ~ismissing(save_to_file)
        save(string(save_to_file, ".jld"), 
             "results", results,
             "f_d_cases", f_d_cases,
             "fd_rate_cases", fd_rate_cases,
             "h_parms", h_parms,
             "Ts", Ts,
             "CN0s", CN0s,
             "state_num", state_num,
             "q_a", qₐ,
             "acquisition_T", acquisition_T,
             "fine_acq_T", fine_acq_T,
             "M", M,
             "dll_b", dll_b,
             "σω", σω,
             "channel", channel,
             "sigtype_name", sigtype.name,
             "nADC", nADC,
             "t_length", t_length,
             "f_s", f_s,
             "doppler_bin_num", doppler_bin_num,
             "p_fa", p_fa,
             "p_d", p_d,
             "f_if", f_if,
             "attempts", attempts,
             "case_names", case_names,
             "use_fd_rate_init", use_fd_rate_init)
    end
    return results
end

function make_plot_for_tracking_cases(file_name; figsize=(8, 3), 
                                      save_to_file=false,
                                      case_names=missing)
    data = load(string(file_name, ".jld"), 
                "results",
                "f_d_cases",
                "fd_rate_cases",
                "h_parms",
                "Ts",
                "CN0s",
                "state_num",
                "q_a")
    results, f_d_cases, fd_rate_cases, h_parms, Ts, CN0s, state_num, qₐ = data
    cmap = get_cmap("viridis")
    fig = figure(figsize=figsize)
    total_subplots = length(Ts)
    cols = 2
    rows = ceil(Int, total_subplots/cols)
    for i in 1:length(Ts)
        T = Ts[i]
        ax = fig.add_subplot(rows, cols, i)
        for j in 1:length(f_d_cases)
            f_d = f_d_cases[j]
            fd_rate = fd_rate_cases[j]
            color = cmap(float(j)/length(f_d_cases))
            if ismissing(case_names)
                label = "Case $(j)"
            else
                label = case_names[j]
            end
            if i == 1
                label = string(label, "; TCXO")
                ax.plot(CN0s, results[i,j,1,:].*(180/π), "-", color=color,
                        label=label)
                ax.plot(CN0s, results[i,j,2,:].*(180/π), ":", color=color)
            elseif i == 2
                label = string(label, "; OCXO")
                ax.plot(CN0s, results[i,j,1,:].*(180/π), "-", color=color)
                ax.plot(CN0s, results[i,j,2,:].*(180/π), ":", color=color,
                        label=label)
            else
                ax.plot(CN0s, results[i,j,1,:].*(180/π), "-", color=color)
                ax.plot(CN0s, results[i,j,2,:].*(180/π), ":", color=color)
            end
        end
        xlabel("C/N₀ (dB⋅Hz)")
        ylabel("Unfiltered σᵩ (degrees)")
        letter = ('a' + (i-1))
        title("($(letter))")
        legend()
    end
    subplots_adjust(wspace=0.35, bottom=0.15, left=0.06, right=0.94, top=0.92)
    if save_to_file
        savefig(string(file_name, ".pdf"), dpi=300)
    end
end


low_fd = 4e3  # GPS max(f_d)
high_fd = 35.2e3  # Starlink max(f_d)

low_fd_rate = -0.730  # GPS min(fd_rate)
high_fd_rate = -475.0  # Starlink min(fd_rate)


f_d_cases     = [     low_fd,       low_fd,      high_fd,     high_fd]
fd_rate_cases = [low_fd_rate, high_fd_rate, high_fd_rate, low_fd_rate]
case_names = [string("Low ", L"f_d", "; Low ", L"\dot{f_d}"),
              string("Low ", L"f_d", "; High ", L"\dot{f_d}"),
              string("High ", L"f_d", "; High ", L"\dot{f_d}"),
              string("High ", L"f_d", "; Low ", L"\dot{f_d}")]
case_num = length(f_d_cases)

tcxo_h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # TCXO
ocxo_h_parms = [2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16] ./ 10e6^2  # OCXO
h_parms = [tcxo_h_parms, ocxo_h_parms]

CN0s = Array(range(20, 50, step=1))

iterations = 10

channel = "I"
sigtype = define_l1ca_code_type()
# Ts = [1e-3 2e-3 5e-3 10e-3]
Ts = [1e-3 10e-3]
f_s = 2.5e6
name = "l1ca"

# channel = "Q"
# sigtype = define_l5_code_type(;channel=channel)
# Ts = [1e-3 20e-3]
# f_s = 25e6
# name = "l5"

t_length = 2
p_fa = 1e-7
qₐ = 0.1
state_num = 2
σω = 1
acquisition_T=20e-3
fine_acq_T=50e-3
M=10
dll_b=0.1


directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name_2state = string(directory, "ch4_compare_tracking_with_diff_cases_$(state_num)-state_q_a=$(qₐ)")

results = compare_tracking_cases(f_d_cases, fd_rate_cases, h_parms, CN0s, Ts,
                                 sigtype, f_s, t_length,
                                 qₐ, state_num, σω, acquisition_T, M,
                                 fine_acq_T, dll_b; 
                                 doppler_bin_num=2, f_if=0, nADC=8,
                                 iterations=iterations, p_fa=1e-7, p_d=0.9,
                                 attempts=3, save_to_file=file_name_2state, 
                                 channel=channel, case_names=case_names) 


t_length = 2
p_fa = 1e-7
qₐ = 0.1
state_num = 2
σω = 1
acquisition_T=20e-3
fine_acq_T=50e-3
M=10
dll_b=0.1


directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name_2state_give_fd_rate = string(directory, "ch4_compare_tracking_with_diff_cases_$(state_num)-state_q_a=$(qₐ)_give_fd_rate")

results = compare_tracking_cases(f_d_cases, fd_rate_cases, h_parms, CN0s, Ts,
                                 sigtype, f_s, t_length,
                                 qₐ, state_num, σω, acquisition_T, M,
                                 fine_acq_T, dll_b; 
                                 doppler_bin_num=2, f_if=0, nADC=8,
                                 iterations=iterations, p_fa=1e-7, p_d=0.9,
                                 attempts=3, 
                                 save_to_file=file_name_2state_give_fd_rate, 
                                 channel=channel, case_names=case_names,
                                 use_fd_rate_init=true) 


t_length = 2
p_fa = 1e-7
qₐ = 0.1
state_num = 3
σω = 1
acquisition_T=20e-3
fine_acq_T=50e-3
M=10
dll_b=0.1

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name_3state_q_a_dont_give_fd_rate = string(directory, "ch4_compare_tracking_with_diff_cases_$(state_num)-state_q_a=$(qₐ)_dont_give_fd_rate")

results = compare_tracking_cases(f_d_cases, fd_rate_cases, h_parms, CN0s, Ts,
                                 sigtype, f_s, t_length,
                                 qₐ, state_num, σω, acquisition_T, M,
                                 fine_acq_T, dll_b; 
                                 doppler_bin_num=2, f_if=0, nADC=8,
                                 iterations=iterations, p_fa=1e-7, p_d=0.9,
                                 attempts=3, 
                                 save_to_file=file_name_3state_q_a_dont_give_fd_rate, 
                                 channel=channel,
                                 case_names=case_names)


t_length = 2
p_fa = 1e-7
qₐ = 0.1
state_num = 3
σω = 1
acquisition_T=20e-3
fine_acq_T=50e-3
M=10
dll_b=0.1

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name_3state_q_a_give_fd_rate = string(directory, "ch4_compare_tracking_with_diff_cases_$(state_num)-state_q_a=$(qₐ)_give_fd_rate")

results = compare_tracking_cases(f_d_cases, fd_rate_cases, h_parms, CN0s, Ts,
                                 sigtype, f_s, t_length,
                                 qₐ, state_num, σω, acquisition_T, M,
                                 fine_acq_T, dll_b; 
                                 doppler_bin_num=2, f_if=0, nADC=8,
                                 iterations=iterations, p_fa=1e-7, p_d=0.9,
                                 attempts=3, 
                                 save_to_file=file_name_3state_q_a_give_fd_rate, 
                                 channel=channel,
                                 case_names=case_names,
                                 use_fd_rate_init=true)


t_length = 2
p_fa = 1e-7
qₐ = 1000
state_num = 3
σω = 1
acquisition_T=20e-3
fine_acq_T=50e-3
M=10
dll_b=0.1

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name_3state = string(directory, "ch4_compare_tracking_with_diff_cases_$(state_num)-state_q_a=$(qₐ)")

results = compare_tracking_cases(f_d_cases, fd_rate_cases, h_parms, CN0s, Ts,
                                 sigtype, f_s, t_length,
                                 qₐ, state_num, σω, acquisition_T, M,
                                 fine_acq_T, dll_b; 
                                 doppler_bin_num=2, f_if=0, nADC=8,
                                 iterations=iterations, p_fa=1e-7, p_d=0.9,
                                 attempts=3, save_to_file=file_name_3state, 
                                 channel=channel,
                                 case_names=case_names)


make_plot_for_tracking_cases(file_name_2state; case_names=case_names,
                             save_to_file=true)
make_plot_for_tracking_cases(file_name_2state_give_fd_rate; case_names=case_names,
                             save_to_file=true)
make_plot_for_tracking_cases(file_name_3state_q_a_dont_give_fd_rate; 
                             case_names=case_names, save_to_file=true)
make_plot_for_tracking_cases(file_name_3state_q_a_give_fd_rate; 
                             case_names=case_names, save_to_file=true)
make_plot_for_tracking_cases(file_name_3state; case_names=case_names,
                             save_to_file=true)