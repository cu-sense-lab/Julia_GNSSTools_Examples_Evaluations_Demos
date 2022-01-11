using GNSSTools
using JLD
using Statistics
using PyPlot
pygui(true)
include("ch4_data_plot_functions.jl")
include("ch4_simulation_function.jl")

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
folder = "figures/l1ca_leo_sat_sim/"

# Iridium
a = Rₑ+780*1000;  # meters
plane_num = 6;
sat_per_plane = 11;
incl = 86;
ΔΩ = 360/plane_num/2;
Δf_per_plane = 360/sat_per_plane/2;  # degrees
eop = get_eop();
CN0 = 45
h_parms = [1.9e-6, 2e-6, 9.7e-8, 9.1e-11, 9.9e-15] ./ 16.3676e6^2  # Rakon IT5300B TCXO
trackresults, data, doppler_curve, doppler_t, sat_elevations, 
sat_azimuth = demo(a, plane_num, sat_per_plane; ΔΩ=ΔΩ, eop=eop, t_length=10*60, 
                   include_phase_noise=true, include_databits=false,
                   h_parms=h_parms, CN0=CN0, showplot=false);

save(string(directory, folder, "ch4_l1ca_leo_sat_data_no_databits.jld"),
     "trackresults", trackresults,
     "data", data,
     "doppler_curve", doppler_curve,
     "doppler_t", doppler_t,
     "sat_elevations", sat_elevations,
     "sat_azimuth", sat_azimuth,
     "a", a,
     "plane_num", plane_num,
     "sat_per_plane", sat_per_plane,
     "incl", incl,
     "d_omega", ΔΩ,
     "Δf_per_plane", Δf_per_plane,
     "h_parms", h_parms,
     "CN0", CN0)

trackresults, data, doppler_curve, doppler_t, sat_elevations, sat_azimuth,
a, plane_num, sat_per_plane, incl, ΔΩ,Δf_per_plane = load(string(directory, folder, "ch4_l1ca_leo_sat_data_no_databits.jld"),
                                                          "trackresults",
                                                          "data",
                                                          "doppler_curve",
                                                          "doppler_t",
                                                          "sat_elevations",
                                                          "sat_azimuth",
                                                          "a",
                                                          "plane_num",
                                                          "sat_per_plane",
                                                          "incl",
                                                          "d_omega",
                                                          "Δf_per_plane")

file_name = string(directory, folder, "ch4_l1ca_leo_sat_example_no_databits.png")
phi0 = π/4
doppler_rate = diff(doppler_curve)
doppler_rate = [doppler_rate[1]; doppler_rate]
n0 = trackresults.data_init_code_chip
cn0 = 45
Tsys = 535
B = 2*l1ca_chipping_rate
A = sqrt(2*GNSSTools.k*Tsys)*10^(cn0/20)
σ = sqrt(GNSSTools.k*B*Tsys)
PS = floor(Int, 5e6*1e-3)*A^2
PN = σ^2
snr = calc_snr(PS, PN)  
plot_example_track_results(trackresults; f_if=data.f_if, 
                           code_length=l1ca_code_length,
                           save_to_file=file_name,
                           bottom=0.075, top=0.95, left=0.125, right=0.875, 
                           hspace=0.75, wspace=0.35, marker_size=1,
                           line_width=0.8, use_ylims=false)#,
                         #   truth_t=t, truth_doppler=doppler, 
                         #   truth_doppler_rate=doppler_rate,
                         #   truth_code_phase=code_phase,
                         #   truth_snr=snr, truth_phi=phi)