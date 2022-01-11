using GNSSTools
using FFTW
using AllanDeviations
using Statistics
using LsqFit
using NNLS
using PyPlot
pygui(true)

function psd_at_f(h_parms, f, v_0)
    h₋₂, h₋₁, h₀, h₁, h₂ = h_parms
    return v_0^2 * (h₋₂*f^(-4) + h₋₁*f^(-3) + h₀*f^(-2) + h₁*f^(-1) + h₂*f^0)
end

function fit_h_parms_to_signal_PSD(psd_dBc, f_c=10e6, 
                                   freqs=[1,10,100,1e3,1e4,1e5])
    N = length(freqs)
    @assert N == length(psd_dBc)
    b = zeros(N-1)
    A = zeros(N-1)
    for i in 1:N-1
        j = N - i
        A[i] = (psd_dBc[j] + psd_dBc[j+1])/2 + 10*log10(freqs[j+1] - freqs[j])
        if i == 1
            b[i] = 2*10^(A[i]/10) / (freqs[j+1] - freqs[j])  # b₀
        elseif i == 2
            b[i] = 2*10^(A[i]/10) / log(freqs[j+1] - freqs[j])  # b₋₁
        else
            b[i] = 2*10^(A[i]/10) / ((2-i)*(freqs[j+1]^(2-i) - freqs[j]^(2-i)))
        end
    end
    # order is h2, h1, h0, h-1, h-2
    h = b ./ f_c.^2
    # reverse order to h-2, h-1, h₀, h₁, h₂
    return reverse(h)
end

function calc_Sy(x, f_s)
    N = length(x)
    X = fft(x)[1:floor(Int, N/2)+1]
    X = (1/(f_s*N)) .* abs2.(X)
    X[2:end-1] .= 2 .* X[2:end-1]
    X_log10 = 10 .* log10.(X)
    frequencies = Array(range(0, f_s/2, length=length(X)))
    return (X, X_log10, frequencies)
end

t_length = 1.  # second
f_s = 5e6  # Hz
v_0 = L1_freq  # Hz
h_parms = h_parms_tcxo[4]
h₋₂, h₋₁, h₀, h₁, h₂ = h_parms

phase_noise = generate_phase_noise(t_length, f_s, v_0, h_parms)

truth = zeros(floor(Int, length(phase_noise)/2)+1)
freqs = Array(range(0, f_s/2, length=length(truth)))
for i in 2:length(truth)
    f = freqs[i]
    truth[i] = v_0^2 * (h₋₂*f^(-4) + h₋₁*f^(-3) + h₀*f^(-2) + h₁*f^(-1) + h₂*f^0)
end

# Take the difference
N = length(phase_noise)
X = fft(phase_noise)[1:floor(Int, N/2)+1]
X = (1/(f_s*N)) .* abs2.(X)
X[2:end-1] .= 2 .* X[2:end-1]

freqs2 = freqs[2:end]
difference = 10 .* log10.(abs.(truth .- X))  # dB
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(freqs, difference, "k")


fig = figure()
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(freqs2, 10 .* log10.(X[2:end]), "b", label="Generated Phase Noise PSD",
         linewidth=5)
ax1.plot(freqs2, 10 .* log10.(truth[2:end]), "r", label="Truth PSD",
         linestyle="--")
ax1.set_xscale("log")
xlabel("Frequency (Hz)")
ylabel("Power")
title("Generated and Truth PSDs")
legend()
ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(freqs2, difference[2:end], "k")
ax2.set_xscale("log")
xlabel("Frequency (Hz)")
ylabel("Power")
title("PSD of the Difference Between the Generated and Truth")
ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(Array(range(0, t_length, length=N)), phase_noise, "k")
xlabel("Time (s)")
ylabel("Amplitude (rad)")
title("Generated Phase Noise in Time Domain")
fig.subplots_adjust(hspace=0.9)


prn = 1
channel = "I"
signal_type = define_l1ca_code_type()
signal = definesignal(signal_type, f_s, t_length; f_d=800, CN0=38.5,
                      receiver_h_parms=h_parms, code_start_idx=100)
generatesignal!(signal)
acqresults, trackresults = process(signal, signal_type, prn, channel;
                                   show_plot=false)
# plotresults(trackresults)
phase_residuals = 2π.*signal.f_d.*trackresults.t .- trackresults.phi
fig = figure()
ax1 = fig.add_subplot(1, 1, 1)
# ax1.plot(trackresults.t, trackresults.phi, "k")
# ax1.plot(trackresults.t, 2π.*signal.f_d.*trackresults.t, "r")
ax1.plot(trackresults.t, phase_residuals, "k")
ax1.plot(signal.t, signal.phase_noise, "b")


################################################################################
f_carrier = 10e6
f_s = 25e6
t_length = 1.
t = Array(0:1/f_s:t_length-1/f_s)
carrier = cis.(2π.*f_carrier.*t)
N = length(t)
N_over_2 = floor(Int, length(carrier)/2)
frequencies = Array(range(0, f_s/2, length=N_over_2+1))
fig = figure()
ax1 = fig.add_subplot(1, 1, 1)
cm_red = get_cmap("Reds")
for i in 1:length(h_parms_tcxo)
    h_parms = h_parms_tcxo[i]
    ϕ_noise = real.(generate_phase_noise(t_length, f_s, f_carrier, h_parms))
    oscillator_noise = cis.(ϕ_noise)
    x = carrier .* oscillator_noise
    X = fft(x)[1:floor(Int, N/2)+1]
    X = (1/(f_s*N)) .* abs2.(X)
	X[2:end-1] .= 2 .* X[2:end-1]
    X_log10 = 10 .* log10.(X)
    max_val = maximum(X_log10)
    color = cm_red(1. * i / length(h_parms_tcxo))
    ax1.plot(frequencies .- f_carrier, X_log10 .- max_val, color=color)
end

cm_blue = get_cmap("Blues")
for i in 1:length(h_parms_ocxo)
    h_parms = h_parms_ocxo[i]
    ϕ_noise = real.(generate_phase_noise(t_length, f_s, f_carrier, h_parms))
    oscillator_noise = cis.(ϕ_noise)
    x = carrier .* oscillator_noise
    X = fft(x)[1:floor(Int, N/2)+1]
    X = (1/(f_s*N)) .* abs2.(X)
	X[2:end-1] .= 2 .* X[2:end-1]
    X_log10 = 10 .* log10.(X)
    max_val = maximum(X_log10)
    color = cm_blue(1. * i / length(h_parms_ocxo))
    ax1.plot(frequencies .- f_carrier, X_log10 .- max_val, color=color)
end

ϕ_noise = real.(generate_phase_noise(t_length, f_s, f_carrier, h_parms_cesium))
oscillator_noise = cis.(ϕ_noise)
x = carrier .* oscillator_noise
X = fft(x)[1:floor(Int, N/2)+1]
X = (1/(f_s*N)) .* abs2.(X)
X[2:end-1] .= 2 .* X[2:end-1]
X_log10 = 10 .* log10.(X)
max_val = maximum(X_log10)
ax1.plot(frequencies .- f_carrier, X_log10 .- max_val, "k--")

xscale("log")
xlim([1, 1e6])
ylim([-200, 0])
xlabel("Freguency (Hz)")
ylabel("dBc/Hz")

# Allan Deviation plot of first four TCXOs
fig = figure()
ax1 = fig.add_subplot(1, 1, 1)
for i in 1:4
    ϕ_noise = real.(generate_phase_noise(t_length, f_s, f_carrier, h_parms_tcxo[i]))
    oscillator_noise = cis.(ϕ_noise)
    x = carrier .* oscillator_noise
    result = allandev(oscillator_noise, f_s)
    ax1.plot(result.tau, result.deviation)
end
τ = Array(range(1e-7, 1, length=1000))
h₋₂, h₋₁, h₀, h₁, h₂ = h_parms_tcxo[4]
σₜ = zeros(length(τ))
for i in 1:length(τ)
    σₜ[i] = sqrt(h₋₂*τ[i]*(2π)^2/6 + h₋₁*2*log(2) + h₀/(2*τ[i]))
end
xscale("log")
yscale("log")

################################################################################
file = "/home/sjbilardi/Documents/SeNSeLab/GNSSTools.jl/test/data/known_oscillator_data/20210327_013800_000005_Aero-North-TCXO_1575.42_25.sc16"
data = loaddata(:sc16, file, 25e6, L1_freq, L1_freq, 1., skip_to=1)
replica = definesignal(define_l1ca_code_type(), 25e6, 1e-3)
for i in 1:32
    fd_est, n0_est, SNR_est = courseacquisition(data, replica, i)
    println("PRN $i\t$SNR_est")
end

# acqresults, trackresults = process(data, define_l1ca_code_type(), 1);

fig = figure()
for i in 1:32
    ax = fig.add_subplot(8, 4, i)
    acqresults, trackresults = process(data, define_l1ca_code_type(), i;
                                       show_plot=false, replica_t_length=1e-3,
                                       fine_acq_method=:carrier, state_num=3,
                                       dynamickf=false)
    ax.plot(trackresults.t, real.(trackresults.ZP))
    ax.plot(trackresults.t, imag.(trackresults.ZP))
    tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    title("PRN $i")
end
subplots_adjust(hspace=0.4, wspace=0.4,  top=0.91, left=0.08, right=0.93, 
                bottom=0.06)

# Track PRN 25
acqresults, trackresults = process(data, define_l1ca_code_type(), 25;
                                   show_plot=false, replica_t_length=1e-3,
                                   fine_acq_method=:carrier)

# Detrend phase and plot PSD
phi = trackresults.phi .- 2π.*median(trackresults.fds).*trackresults.t
# plot_spectrum(phi, 1/replica.t_length)

f_s = 1e3
N = length(phi)
X = fft(phi)[1:floor(Int, N/2)+1]
X = (1/(f_s*N)) .* abs2.(X)
X[2:end-1] .= 2 .* X[2:end-1]
X_phi = deepcopy(X)
frequencies = Array(range(0, f_s/2, length=floor(Int, N/2)+1))
frequencies_phi = deepcopy(frequencies)
v_0 = L1_freq
@. model(f, p) = v_0^2 * (p[1]*f^(-4) + p[2]*f^(-3) + p[3]*f^(-2) + p[4]*f^(-1) + p[5]*f^0)
@. model_truncated(f, p) = v_0^2 * (p[1]*f^(-4) + p[2]*f^(-3) + p[3]*f^(-2))
p₀ = h_parms_tcxo[4]
p₀2 = h_parms_tcxo[1]
# fit = curve_fit(model, frequencies[2:end], X[2:end], p₀)
fit = curve_fit(model_truncated, frequencies[2:end], X_phi[2:end], p₀2)

f_s = 5e6
t_length = 1.
# f_carrier = 10e6
# t_length = 4

# ϕ_noise = real.(generate_phase_noise(t_length, f_s, f_carrier, fit.param))
# t = Array(range(0, t_length, length=length(ϕ_noise)))
# carrier = cis.(2π.*f_carrier.*t)
# oscillator_noise = cis.(ϕ_noise)
# N = length(ϕ_noise)
# x = carrier .* oscillator_noise
# X = fft(x)[1:floor(Int, N/2)+1]
# X = (1/(f_s*N)) .* abs2.(X)
# X[2:end-1] .= 2 .* X[2:end-1]
# X_log10 = 10 .* log10.(X)
# max_val = maximum(X_log10)
frequencies = Array(range(0, f_s/2, length=floor(Int, f_s*t_length)))

truth = zeros(length(frequencies))
truth2 = zeros(length(frequencies))
for i in 2:length(frequencies)
    truth[i] = psd_at_f(fit.param, frequencies[i], L1_freq)
    truth2[i] = psd_at_f(fit.param, frequencies[i], 10e6)
end
truth = 10 .* log10.(truth[2:end])
truth2 = 10 .* log10.(truth2[2:end])

fig = figure()
ax1 = fig.add_subplot(1,1,1)
ax1.plot(frequencies_phi, 10 .* log10.(X_phi), label="Phase measurements from PLL")
ax1.plot(frequencies[2:end], truth, label="Three h-parameter fit")
ax1.plot(frequencies[2:end], truth2, label="Est. of phase noise at 10MHz")
xlim([1, 1e6])
xscale("log")
xlabel("Freguency (Hz)")
ylabel("dB/Hz")
legend()
# ax2 = fig.add_subplot(2,1,2)
# ax2.plot(frequencies .- f_carrier, X_log10 .- max_val, "k-")
# xscale("log")
# xlim([1, 1e6])
# ylim([-200, 0])
# xlabel("Freguency (Hz)")
# ylabel("dBc/Hz")

################################################################################\
# f_s = 36e6
f_carrier1 = 16.3676e6
f_carrier2 = 10e6
f_s = 1*f_carrier1
t_length = 1
N = floor(Int, t_length*f_s)
N_over_2 = floor(Int, N/2)
t = Array(range(0, t_length, length=N))
frequencies = Array(range(0, f_s/2, length=N_over_2+1))
# carrier1 = cis.(2π.*f_carrier1.*t)
# carrier2 = cis.(2π.*f_carrier2.*t)
phase_noise = zeros(Complex{Float64}, N)
# signal = zeros(Complex{Float64}, N)
dBc_freq = [1., 1e1, 1e2, 1e3, 1e4]
X = zeros(Complex{Float64}, N_over_2+1)
x = zeros(Complex{Float64}, N)
rakon_it5300B = Dict("h_parms"=>[1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15],
                     "dBc"=>[-57, -88, -112, -130, -140],
                     "dBc_freq"=>[1., 1e1, 1e2, 1e3, 1e4],
                     "name"=>"TCXO, Rakon IT5300B",
                     "f_0"=>f_carrier1)
agilent_e4424b = Dict("h_parms"=>[2.8e-9, 1.8e-11, 1.2e-10, 1e-11, 1.8e-15],
                      "dBc"=>[-88, -118, -128, -135, -147],
                      "dBc_freq"=>[1., 1e1, 1e2, 1e3, 1e4],
                      "name"=>"OCXO, Agilent E4424B",
                      "f_0"=>10e6)
isotemp_91_1 = Dict("h_parms"=>[2.2e-9, 1.6e-10, 1.6e-10, 7.1e-13, 7.2e-16],
                    "dBc"=>[-89, -120, -140, -151, -154],
                    "dBc_freq"=>[1., 1e1, 1e2, 1e3, 1e4],
                    "name"=>"OCXO, ISOTEMP 91-1",
                    "f_0"=>f_carrier2)
oscilloquartz_8607 = Dict("h_parms"=>[3.9e-14, 8.7e-13, 2.2e-14, 3.2e-13, 6.4e-15],
                          "dBc"=>[-122, -137, -143, -145, -145],
                          "dBc_freq"=>[1., 1e1, 1e2, 1e3, 1e4],
                          "name"=>"OCXO, Oscilloquartz 8607",
                          "f_0"=>f_carrier2)
srs_prs10 = Dict("h_parms"=>[5.3e-11, 4.7e-11, 1.6e-15, 9.6e-14, 9.9e-16],
                 "dBc"=>[-103, -135, -150, -152, -153],
                 "dBc_freq"=>[1., 1e1, 1e2, 1e3, 1e4],
                 "name"=>"Rubidum, SRS PRS10",
                 "f_0"=>f_carrier2)
cesium = Dict("h_parms"=>[1.5e-12, 6.2e-11, 7.6e-23, 1.4e-14, 6.2e-19],
              "dBc"=>[-105, -135, -160, -170, -190],
              "dBc_freq"=>[1., 1e1, 1e2, 1e3, 1e4],
              "name"=>"Cesium, Typical Values",
              "f_0"=>f_carrier2)
mti_250l = Dict("h_parms"=>[1.7e-9, 1e-10, 1.6e-10, 3.5e-13, 1.8e-16],
              "dBc"=>[-90, -120, -140, -155, -160],
              "dBc_freq"=>[1., 1e1, 1e2, 1e3, 1e4],
              "name"=>"OCXO, MTI 250L",#, Low g-sensitivity",
              "f_0"=>f_carrier2)

oscillators = [rakon_it5300B, agilent_e4424b, isotemp_91_1, oscilloquartz_8607, 
               srs_prs10, cesium, mti_250l]
# figsize = (6.5, 4)
colorz = ["k", "b", "g", "r", "y", "m", "c", "gray"]
fig = figure()
ax1 = fig.add_subplot(1,1,1)
for i in 1:length(oscillators)
    color = colorz[i]
    oscillator = oscillators[i]
    if oscillator["name"] == "TCXO, Rakon IT5300B"
        # carrier = carrier1
        f_carrier = f_carrier1
    else
        # carrier = carrier2
        f_carrier = f_carrier2
    end
    h_parms = oscillator["h_parms"] ./ oscillator["f_0"]^2
    generate_phase_noise!(phase_noise, t_length, f_carrier, h_parms)
    ϕ = real.(phase_noise)
    # signal = carrier .* cis.(ϕ)
    # signal .= carrier .* cis.(real.(phase_noise))
    # x = signal
    x .= real.(phase_noise) .+ 0im
    # X, X_log, freqs = calc_Sy(x, f_s)
    fft!(x);
    X .= abs2.(x[1:N_over_2+1]) / (f_s*N);
    X[2:end-1] .= 2 .* X[2:end-1];
    X .= 10 .* log10.(X);
    # X_max = maximum(real.(X));
    # X = X .- X_max;
    # ax1.plot(frequencies[2:end] .- f_carrier, X[2:end], color, 
    #          label=string(oscillator["name"], " from h Parms"))
    ax1.plot(frequencies[2:end], real.(X[2:end]), color, 
             label=string(oscillator["name"], " from h Parms"))
    ax1.scatter(dBc_freq, oscillator["dBc"].+(10*log10(2)), 
                c=color, label=string(oscillator["name"], " Truth"))
end
xscale("log")
xlim([1, 1e6])
xlabel("Frequency (Hz)")
ylabel("dBc/Hz")
legend()


@. function model(x, p)
    f = x[1]
    y = x[2]
    return  y - p[1]*f^(-4) + p[2]*f^(-3) + p[3]*f^(-2) + p[4]*f^(-1) + p[5]*f^0
end

@. model(x, p) = v_0^2 * (p[1]*f^(-4) + p[2]*f^(-3) + p[3]*f^(-2) + p[4]*f^(-1) + p[5]*f^0)