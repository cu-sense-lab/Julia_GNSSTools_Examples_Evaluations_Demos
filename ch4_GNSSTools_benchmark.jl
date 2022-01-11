using GNSSTools
using Statistics
using ProgressMeter
using BenchmarkTools
using JLD
using PyPlot
pygui(true)

function signal_simulation_run(signal)
    definesignal!(signal; new_phase_noise=true, new_thermal_noise=true)
    generatesignal!(signal)
    return signal
end

function replica_simulation_run(signal)
    definereplica!(signal)
    generatereplica!(signal)
    return signal
end

iterations = 5
include_phase_noise = true
include_thermal_noise = true
include_adc = true
include_carrier = true
# sim_t_lengths = [1e-3, 10e-3, 100e-3, 1, 10]
t_lengths = [1e-3, 10e-3, 20e-3, 100e-3, 500e-3, 1]
code_nums = Array(range(1, 5, step=1))
# replica_t_lengths = [1e-3; Array(range(2e-3, 20e-3, step=2e-3))]
f_ss = Array(range(2*l1ca_chipping_rate, 25*l1ca_chipping_rate, 
                   step=2*l1ca_chipping_rate))
signal_btimes = Array{Float64}(undef, iterations)
replica_btimes = Array{Float64}(undef, iterations)
signal_sim_results = Array{Float64}(undef, length(t_lengths), 
                                    length(code_nums), length(f_ss))
replica_sim_results = Array{Float64}(undef, length(t_lengths), 
                                    length(code_nums), length(f_ss))


p = Progress(iterations*length(t_lengths)*length(code_nums)*length(f_ss), 
             1, "Estimating Sim Runtime...")
for i in 1:length(t_lengths)
    for j in 1:length(code_nums)
        code_num = code_nums[j]
        chipping_rates = fill(l1ca_chipping_rate, code_num)
        codes = Array{Dict{Int, Vector{Int}}}(undef, code_num)
        for number_of_codes in 1:code_num
            codes[number_of_codes] = Dict(1 => l1ca_codes[number_of_codes])
        end
        sigtype = definesignaltype(definecodetype(codes, chipping_rates), L1_freq)
        for k in 1:length(f_ss)
            t_length = t_lengths[i]
            f_s = f_ss[k]
            signal = definesignal(sigtype, f_s, t_length; 
                                  include_phase_noise=include_phase_noise, 
                                  include_thermal_noise=include_thermal_noise, 
                                  include_adc=include_adc, 
                                  include_carrier=include_carrier)
            for iteration in 1:iterations
                signal_btimes[iteration] = @elapsed signal_simulation_run(signal)
                replica_btimes[iteration] = @elapsed replica_simulation_run(signal)
                next!(p)
            end
            signal_sim_results[i,j,k] = mean(signal_btimes[2:end])
            replica_sim_results[i,j,k] = mean(replica_btimes[2:end])
        end
    end
end

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name = "ch4_GNSSTools_benchmark_simulation"
file = string(directory, file_name, ".jld")
save(file, 
     "signal_sim_results", signal_sim_results, 
     "replica_sim_results", replica_sim_results,
     "t_lengths", t_lengths,
     "code_nums", code_nums,
     "f_ss", f_ss)


##################### Signal and Replica Simulation times ######################
### Figures for Signal Simulation
directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name = "ch4_GNSSTools_benchmark_simulation"
file = string(directory, file_name, ".jld")
signal_sim_results, replica_sim_results, t_lengths, code_nums, f_ss  = load(
                                               file, 
                                               "signal_sim_results", 
                                               "replica_sim_results",
                                               "t_lengths",
                                               "code_nums",
                                               "f_ss")
fig = figure(figsize=(8.5, 7.5))
cmap = get_cmap("viridis")
ax1 = fig.add_subplot(2, 2, 1)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    ax1.plot(t_lengths, signal_sim_results[:,1,i],
             color=color, label="fₛ = $(round(f_ss[i]*1e-6, digits=1)) MHz")
end
ax1.hlines(y=maximum(signal_sim_results[:,1,:]), xmin=t_lengths[1], 
           xmax=t_lengths[end], linestyles="dotted", 
           color="grey")
xlabel("Simulated Signal Length\n with $(code_nums[1]) Ranging Code (Seconds)")
ylabel("Simulation Runtime (Seconds)")
xlim([t_lengths[1], t_lengths[end]])
ylim([0, maximum(signal_sim_results)])
title("(a)")
ax2= fig.add_subplot(2, 2, 2)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    ax2.plot(t_lengths, signal_sim_results[:,end,i], color=color)
end
ax2.hlines(y=maximum(signal_sim_results[:,1,:]), xmin=t_lengths[1], 
           xmax=t_lengths[end], linestyles="dotted", 
           color="grey", label="Max runtime for 1 ranging code")
# ax2.legend(loc="upper left")
xlabel("Simulated Signal Length\n with $(code_nums[end]) Ranging Codes (Seconds)")
ylabel("Simulation Runtime (Seconds)")
xlim([t_lengths[1], t_lengths[end]])
ylim([0, maximum(signal_sim_results)])
title("(b)")
subplots_adjust(wspace=0.3, bottom=0.1, left=0.08, right=0.92, top=0.93)
### Figures for Replica Simulation
ax3 = fig.add_subplot(2, 2, 3)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    label = string(L"$f_s = $", "$(ceil(Int, f_ss[i]/l1ca_chipping_rate))",
                    L"f_{code}")
    # ax3.plot(t_lengths, replica_sim_results[:,1,i],
    #          color=color, label="fₛ = $(round(f_ss[i]/l1ca_chipping_rate, digits=0)) MHz")
    ax3.plot(t_lengths, replica_sim_results[:,1,i],
             color=color, label=label)
end
ax3.hlines(y=maximum(replica_sim_results[:,1,:]), xmin=t_lengths[1], 
           xmax=t_lengths[end], linestyles="dotted", 
           color="grey")
ax3.legend(loc="upper left", ncol=2)
xlabel("Generated Replica Length\n with $(code_nums[1]) Ranging Code (Seconds)")
ylabel("Simulation Runtime (Seconds)")
xlim([t_lengths[1], t_lengths[end]])
ylim([0, maximum(replica_sim_results)])
title("(c)")
ax4= fig.add_subplot(2, 2, 4)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    ax4.plot(t_lengths, replica_sim_results[:,end,i], color=color)
end
ax4.hlines(y=maximum(replica_sim_results[:,1,:]), xmin=t_lengths[1], 
           xmax=t_lengths[end], linestyles="dotted", 
           color="grey", label="Max runtime for 1 ranging code")
ax4.legend(loc="upper left")
xlabel("Generated Replica Length\n with $(code_nums[end]) Ranging Codes (Seconds)")
ylabel("Simulation Runtime (Seconds)")
xlim([t_lengths[1], t_lengths[end]])
ylim([0, maximum(replica_sim_results)])
title("(d)")
subplots_adjust(wspace=0.3, hspace=0.45, bottom=0.12, 
                left=0.08, right=0.92, top=0.93)
savefig(string(directory, "figures/", file_name, ".svg"), dpi=300)


############################## Course Acquisition ############################## 
function course_acquisition_run(signal, replica, M)
    courseacquisition(signal, replica, 1; M=M)
end

iterations = 3
# sim_t_lengths = [1e-3, 10e-3, 100e-3, 1, 10]
t_lengths = [1e-3, 2e-3, 5e-3, 10e-3,  20e-3, 30e-3, 40e-3, 50e-3]
Ms = floor.(Int, t_lengths./1e-3)
code_nums = Array(range(1, 5, step=1))
# replica_t_lengths = [1e-3; Array(range(2e-3, 20e-3, step=2e-3))]
f_ss = Array(range(2*l1ca_chipping_rate, 25*l1ca_chipping_rate, 
                   step=2*l1ca_chipping_rate))
acquisition_coherent_btimes = Array{Float64}(undef, iterations)
acquisition_noncoherent_btimes = Array{Float64}(undef, iterations)
acquisition_coherent_results = Array{Float64}(undef, length(t_lengths), 
                                              length(code_nums), length(f_ss))
acquisition_noncoherent_results = Array{Float64}(undef, length(t_lengths), 
                                                 length(code_nums), length(f_ss))

p = Progress(iterations*length(t_lengths)*length(code_nums)*length(f_ss), 
             1, "Estimating Acquisition Runtime...")
for i in length(t_lengths):-1:1
    for j in length(code_nums):-1:1
        code_num = code_nums[j]
        chipping_rates = fill(l1ca_chipping_rate, code_num)
        codes = Array{Dict{Int, Vector{Int}}}(undef, code_num)
        for number_of_codes in 1:code_num
            codes[number_of_codes] = Dict(1 => l1ca_codes[number_of_codes])
        end
        sigtype = definesignaltype(definecodetype(codes, chipping_rates), L1_freq)
        for k in length(f_ss):-1:1
            t_length = t_lengths[i]
            M = Ms[i]
            f_s = f_ss[k]
            signal = definesignal(sigtype, f_s, t_length; 
                                  include_phase_noise=false, 
                                  include_thermal_noise=false, 
                                  include_adc=false, 
                                  include_carrier=false,
                                  allocate_noise_vectors=false)
            generatesignal!(signal)
            replica_coherent = definereplica(sigtype, f_s, t_length)
            replica_noncoherent = definereplica(sigtype, f_s, 1e-3)
            for iteration in 1:iterations
                acquisition_coherent_btimes[iteration] = @elapsed course_acquisition_run(
                                                             signal,
                                                             replica_coherent,
                                                             1)
                acquisition_noncoherent_btimes[iteration] = @elapsed course_acquisition_run(
                                                             signal,
                                                             replica_noncoherent,
                                                             M)
                next!(p)
            end
            acquisition_coherent_results[i,j,k] = mean(acquisition_coherent_btimes[2:end])
            acquisition_noncoherent_results[i,j,k] = mean(acquisition_noncoherent_btimes[2:end])
        end
    end
end

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name = "ch4_GNSSTools_benchmark_acquisition"
file = string(directory, file_name, ".jld")
save(file, 
     "acquisition_coherent_results", acquisition_coherent_results, 
     "acquisition_noncoherent_results", acquisition_noncoherent_results,
     "t_lengths", t_lengths,
     "Ms", Ms,
     "code_nums", code_nums,
     "f_ss", f_ss)

################## Plots
directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name = "ch4_GNSSTools_benchmark_acquisition"
file = string(directory, file_name, ".jld")
results = load(file, 
               "acquisition_coherent_results",
               "acquisition_noncoherent_results",
               "t_lengths",
               "Ms",
               "code_nums",
               "f_ss")
acquisition_coherent_results, acquisition_noncoherent_results,
t_lengths, Ms, code_nums, f_ss = results
t_lengths = t_lengths .* 1e3

fig = figure(figsize=(8.5, 7.5))
cmap = get_cmap("viridis")
ax1 = fig.add_subplot(2, 2, 1)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    label = string(L"$f_s = $", "$(ceil(Int, f_ss[i]/l1ca_chipping_rate))",
                    L"f_{code}")
    ax1.plot(t_lengths, acquisition_coherent_results[:,1,i],
             color=color, label=label)
end
ax1.hlines(y=maximum(acquisition_coherent_results[:,1,:]), xmin=t_lengths[1], 
           xmax=t_lengths[end], linestyles="dotted", 
           color="grey")
xlabel("Integration Time\n with $(code_nums[1]) Ranging Code (ms)")
ylabel("Simulation Runtime (Seconds)")
xlim([t_lengths[1], t_lengths[end]])
ylim([0, maximum(acquisition_coherent_results)])
title("(a)")
ax1.legend(loc="upper left", ncol=2)
ax2= fig.add_subplot(2, 2, 2)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    ax2.plot(t_lengths, acquisition_coherent_results[:,end,i], color=color)
end
ax2.hlines(y=maximum(acquisition_coherent_results[:,1,:]), xmin=t_lengths[1], 
           xmax=t_lengths[end], linestyles="dotted", 
           color="grey", label="Max runtime for 1 ranging code")
# ax2.legend(loc="upper left")
xlabel("Integration Time\n with $(code_nums[end]) Ranging Codes (ms)")
ylabel("Simulation Runtime (Seconds)")
xlim([t_lengths[1], t_lengths[end]])
ylim([0, maximum(acquisition_coherent_results)])
title("(b)")
### Figures for Replica Simulation
ax3 = fig.add_subplot(2, 2, 3)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    label = string(L"$f_s = $", "$(ceil(Int, f_ss[i]/l1ca_chipping_rate))",
                    L"f_{code}")
    # ax3.plot(t_lengths, replica_sim_results[:,1,i],
    #          color=color, label="fₛ = $(round(f_ss[i]/l1ca_chipping_rate, digits=0)) MHz")
    ax3.plot(Ms, acquisition_noncoherent_results[:,1,i],
             color=color, label=label)
end
ax3.hlines(y=maximum(acquisition_noncoherent_results[:,1,:]), xmin=Ms[1], 
           xmax=Ms[end], linestyles="dotted", 
           color="grey")
# ax3.legend(loc="upper left", ncol=2)
xlabel("Number of 1ms Integrations\n with $(code_nums[1]) Ranging Code")
ylabel("Simulation Runtime (Seconds)")
xlim([Ms[1], Ms[end]])
ylim([0, maximum(acquisition_noncoherent_results)])
title("(c)")
ax4= fig.add_subplot(2, 2, 4)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    ax4.plot(Ms, acquisition_noncoherent_results[:,end,i], color=color)
end
ax4.hlines(y=maximum(acquisition_noncoherent_results[:,1,:]), xmin=Ms[1], 
           xmax=Ms[end], linestyles="dotted", 
           color="grey", label="Max runtime for 1 ranging code")
ax4.legend(loc="upper left")
xlabel("Number of 1ms Integrations\n with $(code_nums[end]) Ranging Codes")
ylabel("Simulation Runtime (Seconds)")
xlim([Ms[1], Ms[end]])
ylim([0, maximum(acquisition_noncoherent_results)])
title("(d)")
subplots_adjust(wspace=0.3, hspace=0.45, bottom=0.12, 
                left=0.08, right=0.92, top=0.93)
savefig(string(directory, "figures/", file_name, ".svg"), dpi=300)

############################### Signal Tracking ############################### 
function tracking_run(signal, replica)
    R = [1.]
    P = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    prn = 1
    ϕ_init = 0.
    fd_init = 0.
    n0_idx_init = 1.
    trackprn(signal, replica, prn, ϕ_init, fd_init, n0_idx_init, P, R)
end

iterations = 3
# sim_t_lengths = [1e-3, 10e-3, 100e-3, 1, 10]
t_lengths = [0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
Ts = [1e-3, 2e-3, 5e-3, 10e-3, 20e-3]
code_nums = Array(range(1, 5, step=1))
# replica_t_lengths = [1e-3; Array(range(2e-3, 20e-3, step=2e-3))]
f_ss = Array(range(2*l1ca_chipping_rate, 25*l1ca_chipping_rate, 
                   step=2*l1ca_chipping_rate))
tracking_btimes = Array{Float64}(undef, iterations)
tracking_results = Array{Float64}(undef, length(t_lengths), 
                                  length(code_nums), length(f_ss),
                                  length(Ts))

p = Progress(iterations*length(t_lengths)*length(code_nums)*length(f_ss)*length(Ts), 
             1, "Estimating Tracking Runtime...")
for i in length(t_lengths):-1:1
    for j in length(code_nums):-1:1
        code_num = code_nums[j]
        chipping_rates = fill(l1ca_chipping_rate, code_num)
        codes = Array{Dict{Int, Vector{Int}}}(undef, code_num)
        for number_of_codes in 1:code_num
            codes[number_of_codes] = Dict(1 => l1ca_codes[number_of_codes])
        end
        sigtype = definesignaltype(definecodetype(codes, chipping_rates), L1_freq)
        for k in length(f_ss):-1:1
            t_length = t_lengths[i]
            M = Ms[i]
            f_s = f_ss[k]
            signal = definesignal(sigtype, f_s, t_length; 
                                  include_phase_noise=false, 
                                  include_thermal_noise=false, 
                                  include_adc=false, 
                                  include_carrier=false,
                                  allocate_noise_vectors=false)
            generatesignal!(signal)
            for g in 1:length(Ts)
                T = Ts[g]
                replica = definereplica(sigtype, f_s, T)
                for iteration in 1:iterations
                    tracking_btimes[iteration] = @elapsed tracking_run(
                                                          signal, 
                                                          replica)
                    next!(p)
                end
                tracking_results[i,j,k,g] = mean(tracking_btimes[2:end])
            end
        end
    end
end

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name = "ch4_GNSSTools_benchmark_tracking"
file = string(directory, file_name, ".jld")
save(file, 
     "tracking_results", tracking_results, 
     "t_lengths", t_lengths,
     "Ms", Ms,
     "code_nums", code_nums,
     "f_ss", f_ss,
     "Ts", Ts)

################## Plots

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
file_name = "ch4_GNSSTools_benchmark_tracking"
file = string(directory, file_name, ".jld")
tracking_results, t_lengths, Ms, code_nums, f_ss, Ts = load(file, 
                                                            "tracking_results",
                                                            "t_lengths",
                                                            "Ms",
                                                            "code_nums",
                                                            "f_ss",
                                                            "Ts")
Ts = Ts .* 1e3


fig = figure(figsize=(8.5, 4))
cmap = get_cmap("viridis")
ax1 = fig.add_subplot(1, 2, 1)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    label = string(L"$f_s = $", "$(ceil(Int, f_ss[i]/l1ca_chipping_rate))",
                    L"f_{code}")
    ax1.plot(Ts, tracking_results[end,1,i,:],
             color=color, label=label)
end
ax1.hlines(y=maximum(tracking_results[end,1,:,:]), xmin=Ts[1], 
           xmax=Ts[end], linestyles="dotted", 
           color="grey")
xlabel("Coherent Integration Time\n with $(code_nums[1]) Ranging Code (ms)")
ylabel("Simulation Runtime (Seconds)")
xlim([Ts[1], Ts[end]])
ylim([0, maximum(tracking_results[end,:,:,:])])
ax1.legend(loc="upper left", ncol=2)
title("(a)")
ax2= fig.add_subplot(1, 2, 2)
for i in 1:length(f_ss)
    color = cmap(float(i)/length(f_ss))
    ax2.plot(Ts, tracking_results[end,end,i,:], color=color)
end
ax2.hlines(y=maximum(tracking_results[end,1,:,:]), xmin=Ts[1], 
           xmax=Ts[end], linestyles="dotted", 
           color="grey", label="Max runtime for 1 ranging code")
# ax2.legend(loc="upper left")
xlabel("Coherent Integration Time\n with $(code_nums[end]) Ranging Codes (ms)")
ylabel("Simulation Runtime (Seconds)")
xlim([Ts[1], Ts[end]])
ylim([0, maximum(tracking_results[end,:,:,:])])
ax2.legend(loc="upper left")
title("(b)")
subplots_adjust(wspace=0.3, hspace=0.55, bottom=0.15, 
                left=0.08, right=0.92, top=0.93)
savefig(string(directory, "figures/", file_name, ".svg"), dpi=300)