using GNSSTools
using FFTW
using AllanDeviations
using Statistics
# using LsqFit
# using NNLS
using PyPlot
pygui(true)

# function psd_at_f(h_parms, f, v_0)
#     h₋₂, h₋₁, h₀, h₁, h₂ = h_parms
#     return v_0^2 * (h₋₂*f^(-4) + h₋₁*f^(-3) + h₀*f^(-2) + h₁*f^(-1) + h₂*f^0)
# end

# function fit_h_parms_to_signal_PSD(psd_dBc, f_c=10e6, 
#                                    freqs=[1,10,100,1e3,1e4,1e5])
#     N = length(freqs)
#     @assert N == length(psd_dBc)
#     b = zeros(N-1)
#     A = zeros(N-1)
#     for i in 1:N-1
#         j = N - i
#         A[i] = (psd_dBc[j] + psd_dBc[j+1])/2 + 10*log10(freqs[j+1] - freqs[j])
#         if i == 1
#             b[i] = 2*10^(A[i]/10) / (freqs[j+1] - freqs[j])  # b₀
#         elseif i == 2
#             b[i] = 2*10^(A[i]/10) / log(freqs[j+1] - freqs[j])  # b₋₁
#         else
#             b[i] = 2*10^(A[i]/10) / ((2-i)*(freqs[j+1]^(2-i) - freqs[j]^(2-i)))
#         end
#     end
#     # order is h2, h1, h0, h-1, h-2
#     h = b ./ f_c.^2
#     # reverse order to h-2, h-1, h₀, h₁, h₂
#     return reverse(h)
# end

function calc_Sy(x, f_s)
    N = length(x)
    X = fft(x)[1:floor(Int, N/2)+1]
    X = (1/(f_s*N)) .* abs2.(X)
    X[2:end-1] .= 2 .* X[2:end-1]
    X_log10 = 10 .* log10.(X)
    frequencies = Array(range(0, f_s/2, length=length(X)))
    return (X, X_log10, frequencies)
end

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
fig = figure(figsize=(6.5,4))
ax1 = fig.add_subplot(1,1,1)
legend_markers = []
marker_names = []
for i in 1:length(oscillators)
# for i in 1:1
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
    # ϕ = real.(phase_noise)
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
    ax1.plot(frequencies[2:end], real.(X[2:end]), color)#), 
            #  label=oscillator["name"])
    # ax1.plot(frequencies[2:end], real.(X[2:end]), color, 
            #  label=string(oscillator["name"], " from h Parms"))
    ax1.scatter(dBc_freq, oscillator["dBc"].+(10*log10(2)), 
                c=color)#, label=string(oscillator["name"], " Truth"))
    legend_marker = plt.Line2D((0,1),(0,0), color=color, marker="o", linestyle="-")
    push!(legend_markers, legend_marker)
    push!(marker_names, oscillator["name"])
end
xscale("log")
xlim([1, 1e6])
xlabel("Frequency (Hz)")
ylabel("dBc/Hz")
subplots_adjust(top=0.96, left=0.12, right=0.88, bottom=0.13)
ax1.legend(legend_markers, marker_names)

# savefig("figures/ch2_oscillator_psds.pdf", dpi=300)