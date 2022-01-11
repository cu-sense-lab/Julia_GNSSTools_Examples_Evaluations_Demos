using GNSSTools
using Statistics
using ProgressMeter
using JLD
using PyPlot
include("ch4_constellation_hist_structs.jl")
pygui(true)

function get_median_doppler(dopplers, doppler_rates, max_elevations, 
                            elev_range; nbins=100)
    elev_min = minimum(elev_range)
    elev_max = maximum(elev_range)
    idxs = (max_elevations .>= elev_min) .== (max_elevations .< elev_max)
    doppler_elev = dopplers[idxs]
    doppler_rate_elev = doppler_rates[idxs]
    x = get_histogram((doppler_rate_elev, doppler_elev); nbins=nbins)
    doppler_elev_mean = x.edges[2][1:end-1] .+ (x.edges[2][2] - x.edges[2][1])/2
    doppler_rate_mean = Array{Float64}(undef, length(doppler_elev_mean))
    for i in 1:size(x.weights)[2]
        # y = x.weights[:,i]
        # edges = x.edges[1][1:end-1]
        # idxs = y .> 0
        # if any(idxs)
        #     doppler_rate_mean[i] = median(edges[y .> 0])
        # else
        #     doppler_rate_mean[i] = NaN
        # end
        doppler_rate_mean[i] = sum(x.weights[:,i].*x.edges[1][1:end-1])/sum(x.weights[:,i])
    end
    return (x, doppler_elev, doppler_rate_elev, 
            doppler_elev_mean, doppler_rate_mean)
end

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
# directory = ""
doppler_file_name = string(directory, "ch1_doppler_curves.jld")
doppler_raw_gps, doppler_rate_raw_gps, 
elevations_gps, sig_freq, max_elevations_gps = load(
                                                       doppler_file_name, 
                                                       "doppler_raw_gps",
                                                       "doppler_rate_raw_gps",
                                                       "elevations_gps",
                                                       "freq",
                                                       "max_elevations_gps")
doppler_raw_iridium, doppler_rate_raw_iridium, 
elevations_iridium, max_elevations_iridium = load(
                                                     doppler_file_name, 
                                                     "doppler_raw_iridium",
                                                     "doppler_rate_raw_iridium",
                                                     "elevations_iridium",
                                                     "max_elevations_iridium")
doppler_raw_starlink, doppler_rate_raw_starlink, 
elevations_starlink, max_elevations_starlink = load(
                                                       doppler_file_name, 
                                                       "doppler_raw_starlink",
                                                       "doppler_rate_raw_starlink",
                                                       "elevations_starlink",
                                                       "max_elevations_starlink")
doppler_raw_oneweb, doppler_rate_raw_oneweb,
elevations_oneweb, max_elevations_oneweb = load(doppler_file_name, 
                         "doppler_raw_oneweb",
                         "doppler_rate_raw_oneweb",
                         "elevations_oneweb",
                         "max_elevations_oneweb")

gps_hist = initiate_hist(doppler_raw_gps, doppler_rate_raw_gps, 
                         elevations_gps, max_elevations_gps; 
                         nbins=1000, elev_bins=100)
iridium_hist = initiate_hist(doppler_raw_iridium, doppler_rate_raw_iridium,
                             elevations_iridium, max_elevations_iridium; 
                             nbins=1000)
starlink_hist = initiate_hist(doppler_raw_starlink, doppler_rate_raw_starlink,
                              elevations_starlink, max_elevations_starlink; 
                              nbins=1000)
oneweb_hist = initiate_hist(doppler_raw_oneweb, doppler_rate_raw_oneweb,
                            elevations_oneweb, max_elevations_oneweb; 
                            nbins=1000)


# elevation_bins = Array(range(1, 91, step=2))
# gps_doppler_max, gps_doppler_max = max_doppler_curve_per_elev(
#                                                     gps_hist, elevation_bins)
# iridium_doppler_max, iridium_doppler_rate_max = max_doppler_curve_per_elev(
#                                                     iridium_hist, elevation_bins)
# starlink_doppler_max, starlink_doppler_rate_max = max_doppler_curve_per_elev(
#                                                     starlink_hist, elevation_bins)
# oneweb_doppler_max, oneweb_doppler_rate_max = max_doppler_curve_per_elev(
#                                                     oneweb_hist, elevation_bins)
                                                    

elev_ranges = [[10, 11],
               [15, 16],
               [30, 31],
               [45, 46],
               [60, 61],
               [90, 89]]
nbins = 100
cmap = "gray_r"
fig = figure(figsize=(7.5, 5))
ax1 = fig.add_subplot(2,2,1)
hist2D(gps_hist.dopplers./1000, gps_hist.doppler_rates, bins=nbins, 
          cmap=cmap, density=true, vmax=1)
# for i in 1:length(elev_ranges)
#     elev_min = elev_ranges[i][1]
#     elev_max = elev_ranges[i][2]
#     x, d, dd, d_mean, dd_mean = get_median_doppler(gps_hist.dopplers, 
#                                                    gps_hist.doppler_rates,
#                                                    gps_hist.max_elevations, 
#                                                    [elev_min, elev_max]; 
#                                                    nbins=(500, 100))
#     ax1.plot(d_mean./1000, dd_mean, label="$(elev_min)ᵒ",
#              linewidth=1)
# end
xlabel("Doppler (kHz)")
ylabel("Doppler Rate (Hz/s)")
# legend()
title("(a)")

ax2 = fig.add_subplot(2,2,2, aspect="auto")
ax2.hist2d(iridium_hist.dopplers./1000, iridium_hist.doppler_rates, bins=nbins, 
            cmap=cmap, density=true, vmax=0.0005)
for i in 1:length(elev_ranges)
    elev_min = elev_ranges[i][1]
    elev_max = elev_ranges[i][2]
    x, d, dd, d_mean, dd_mean = get_median_doppler(iridium_hist.dopplers, 
                                                   iridium_hist.doppler_rates,
                                                   iridium_hist.max_elevations, 
                                                   [elev_min, elev_max]; 
                                                   nbins=(500, 100))
    ax2.plot(d_mean./1000, dd_mean, label="$(elev_min)°",
             linewidth=1)
end
xlabel("Doppler (kHz)")
ylabel("Doppler Rate (Hz/s)")
legend(bbox_to_anchor=(1, 1.04), loc="upper left", ncol=1)
title("(b)")

ax3 = fig.add_subplot(2,2,3, aspect="auto")
ax3.hist2d(starlink_hist.dopplers./1000, starlink_hist.doppler_rates, bins=nbins, 
           cmap=cmap, density=true, vmax=0.0005)
for i in 1:length(elev_ranges)
    elev_min = elev_ranges[i][1]
    elev_max = elev_ranges[i][2]
    x, d, dd, d_mean, dd_mean = get_median_doppler(starlink_hist.dopplers, 
                                                   starlink_hist.doppler_rates,
                                                   starlink_hist.max_elevations, 
                                                   [elev_min, elev_max]; 
                                                   nbins=(500, 100))
    ax3.plot(d_mean./1000, dd_mean,
             linewidth=1)
end
xlabel("Doppler (kHz)")
ylabel("Doppler Rate (Hz/s)")
# legend()
title("(c)")

ax4 = fig.add_subplot(2,2,4, aspect="auto")
ax4.hist2d(oneweb_hist.dopplers./1000, oneweb_hist.doppler_rates, bins=nbins, 
           cmap=cmap, density=true, vmax=0.0005)
for i in 1:length(elev_ranges)
    elev_min = elev_ranges[i][1]
    elev_max = elev_ranges[i][2]
    x, d, dd, d_mean, dd_mean = get_median_doppler(oneweb_hist.dopplers, 
                                                   oneweb_hist.doppler_rates,
                                                   oneweb_hist.max_elevations, 
                                                   [elev_min, elev_max]; 
                                                   nbins=(500, 100))
    ax4.plot(d_mean./1000, dd_mean,
             linewidth=1)
end
xlabel("Doppler (kHz)")
ylabel("Doppler Rate (Hz/s)")
# legend()
title("(d)")

subplots_adjust(wspace=0.4, hspace=0.5, 
                bottom=0.1, left=0.12, right=0.88, top=0.93)

savefig(string(directory, "figures/ch4_joint_doppler-doppler_rate_distributions.pdf"), dpi=300)

fig = figure()
ax = fig.add_subplot(1,1,1)
elev_ranges = [[10, 10.5],
               [15, 15.1],
               [30, 30.1],
               [45, 45.1],
               [60, 60.1],
               [90, 89.9]]
ax.hist2d(doppler_raw_iridium./1000, doppler_rate_raw_iridium, bins=100, 
          cmap="gray_r", density=true, vmax=0.0005)
for i in 1:length(elev_ranges)
    elev_min = elev_ranges[i][1]
    elev_max = elev_ranges[i][2]
    x, d, dd, d_mean, dd_mean = get_median_doppler(doppler_raw_iridium, 
                                                   doppler_rate_raw_iridium,
                                                   max_elevations_iridium, 
                                                   [elev_min, elev_max]; 
                                                   nbins=(500, 100))
    ax.plot(d_mean./1000, dd_mean, label="$(elev_min)ᵒ",
            linewidth=1)
end
xlabel("Doppler (kHz)")
ylabel("Doppler Rate (Hz/s)")
legend()


hist_data = starlink_hist
fig = figure(figsize=(7.5, 4))
nbins = 100
max_elev_lim = 10
max_elev_idxs = hist_data.max_elevations .>= max_elev_lim
ax1 = fig.add_subplot(2, 3, 1)
idxs = (hist_data.elevations .>= 10) .== (hist_data.elevations .< 15) .== max_elev_idxs
ax1.hist(hist_data.dopplers[idxs]./1000,
         color="k", bins=nbins, density=true);
xlabel("Doppler (kHz)")
ylabel("Probability")
xlim([minimum(hist_data.dopplers), maximum(hist_data.dopplers)]./1000)
ax1.ticklabel_format(axis="y", style="sci", scilimits=[-1,1])
title("(a)")

ax2 = fig.add_subplot(2, 3, 2)
idxs = (hist_data.elevations .>= 60) .== (hist_data.elevations .< 65) .== max_elev_idxs
ax2.hist(hist_data.dopplers[idxs]./1000,
         color="k", bins=nbins, density=true);
xlabel("Doppler (kHz)")
ylabel("Probability")
xlim([minimum(hist_data.dopplers), maximum(hist_data.dopplers)]./1000)
ax2.ticklabel_format(axis="y", style="sci", scilimits=[-1,1])
title("(b)")

ax3 = fig.add_subplot(2, 3, 3)
idxs = (hist_data.elevations .>= 85) .== (hist_data.elevations .< 90) .== max_elev_idxs
ax3.hist(hist_data.dopplers[idxs]./1000,
color="k", bins=nbins, density=true);
xlabel("Doppler (kHz)")
ylabel("Probability")
xlim([minimum(hist_data.dopplers), maximum(hist_data.dopplers)]./1000)
ax3.ticklabel_format(axis="y", style="sci", scilimits=[-1,1])
title("(c)")

ax1 = fig.add_subplot(2, 3, 4)
idxs = (hist_data.elevations .>= 10) .== (hist_data.elevations .< 15) .== max_elev_idxs
ax1.hist(hist_data.doppler_rates[idxs],
         color="k", bins=nbins, density=true);
xlabel("Doppler Rate (Hz/s)")
ylabel("Probability")
xlim([minimum(hist_data.doppler_rates), maximum(hist_data.doppler_rates)])
ax1.ticklabel_format(axis="y", style="sci", scilimits=[-1,1])
title("(d)")

ax2 = fig.add_subplot(2, 3, 5)
idxs = (hist_data.elevations .>= 60) .== (hist_data.elevations .< 65) .== max_elev_idxs
ax2.hist(hist_data.doppler_rates[idxs],
         color="k", bins=nbins, density=true);
xlabel("Doppler Rate (Hz/s)")
ylabel("Probability")
xlim([minimum(hist_data.doppler_rates), maximum(hist_data.doppler_rates)])
ax2.ticklabel_format(axis="y", style="sci", scilimits=[-1,1])
title("(e)")

ax3 = fig.add_subplot(2, 3, 6)
idxs = (hist_data.elevations .>= 85) .== (hist_data.elevations .< 90) .== max_elev_idxs
ax3.hist(hist_data.doppler_rates[idxs],
color="k", bins=nbins, density=true);
xlabel("Doppler Rate (Hz/s)")
ylabel("Probability")
xlim([minimum(hist_data.doppler_rates), maximum(hist_data.doppler_rates)])
ax3.ticklabel_format(axis="y", style="sci", scilimits=[-1,1])
title("(f)")

subplots_adjust(wspace=0.36, hspace=0.6, bottom=0.115, 
                left=0.07, right=0.93, top=0.94)
savefig(string(directory, "figures/ch4_starlink_hist_at_different_elevations_above_$(max_elev_lim)_degrees_elevation.pdf"), dpi=300)


######### Skyplots for diff passes with max elevations
# Iridium
δt = 1  # seconds
total_time = 1*24*60*60  # seconds
t_range = Array(0:δt:total_time);
min_elevation = 0  # degrees
user_lla = (40.01, -105.2437, 1655)
a = Rₑ+780*1000;  # meters
plane_num = 6;
sat_per_plane = 11;
incl = 86;
ΔΩ = 360/plane_num/2;
Δf_per_plane = 360/sat_per_plane/2;  # degrees
orbital_period = 2π*sqrt(a^3 / GM)  # seconds
iridium_constellation= define_constellation(a, plane_num, 
                                            sat_per_plane, incl, 
                                            t_range; obs_lla=user_lla, 
                                            show_plot=false,
                                            print_steps=true, 
                                            ΔΩ=ΔΩ,
                                            Δf_per_plane=Δf_per_plane);



elev_azimuths = []
elev_elevations = []
elev_max_elevations = []
for satellite in iridium_constellation.satellites
    ts = Array{Float64}(undef, length(satellite.t))
    ts[1] = 0.0
    azimuths = Array{Float64}(undef, length(satellite.t))
    elevations = Array{Float64}(undef, length(satellite.t))
    max_elevations = Array{Float64}(undef, length(satellite.t))
    k_start = 1
    k = 1
    pass = 1
    pass_num = 3
    for j in 1:length(satellite.t)
        sat_range, azimuth, elevation = calcelevation(satellite.r_ecef[j,:], 
                                                        user_lla)
        if k > 1
            if elevation > min_elevation
                ts[k] = satellite.t[j]
                azimuths[k] = azimuth
                elevations[k] = elevation
                Δt = abs(ts[k] - ts[k-1])*(60*60*24)
                if Δt > 2*δt
                    # if pass == pass_num
                        max_elevations[k_start:k-1] .= maximum(elevations[k_start:k-1])
                        k_start = k
                    # else
                    #     max_elevations[k_start:j-1] .= -1
                    #     k_start = k
                    # end
                    pass += 1
                end
                k += 1
            end
        else
            if elevation > min_elevation
                ts[k] = satellite.t[j]
                azimuths[k] = azimuth
                elevations[k] = elevation
                k += 1
            end
        end
    end
    if k_start < k
        max_elevations[k_start:k-1] .= -1
    end
    push!(elev_azimuths, azimuths[1:k-1])
    push!(elev_elevations, elevations[1:k-1])
    push!(elev_max_elevations, max_elevations[1:k-1])
end

δ_deg = 1
elev_ranges = [[10, 10+δ_deg],
               [15, 15+δ_deg],
               [30, 30+δ_deg],
               [45, 45+δ_deg],
               [60, 60+δ_deg],
               [90, 90-δ_deg]]
fig = figure(figsize=(7.5, 4))
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
          "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
for i in 1:length(elev_ranges)
    ax = fig.add_subplot(2, 3, i, polar=true)
    color = colors[i]
    elev_min = minimum(elev_ranges[i])
    elev_max = maximum(elev_ranges[i])
    for j in 1:length(elev_max_elevations)
        max_els = elev_max_elevations[j]
        idxs = (max_els .>= elev_min) .== (max_els .<= elev_max)
        if any(idxs)
            ax.plot(elev_azimuths[j][idxs].*(π/180), 
                       90 .- elev_elevations[j][idxs], 
                       linestyle="-", linewidth=1, color=color)
            break
        end
    end
    ax.set_theta_zero_location("N")
    # Reverse direction of azimuth to clockwise
    ax.set_theta_direction(-1)
    # Reverse y-axis so that 90 degrees elevation is at center
    ax.set_rlim(0, 90, 1)
    ax.set_yticks(range(0, 91, step=30))
    fontsize = 8
    ax.set_yticklabels(reverse(ax.get_yticks()), fontsize=fontsize-2)
    ax.set_xticklabels(["N", "", "E", "", "S", "", "W", ""], fontsize=fontsize)
    # ylabel("Elevation (deg)")
    letter = 'a' + (i-1)
    title("($letter)")
    # title("$(elev_ranges[i][1])°")
end
subplots_adjust(wspace=0.7, hspace=0.8, bottom=0.07, 
                left=0.02, right=0.98, top=0.88)
savefig(string(directory, "figures/ch4_diff_max_elevation_passes.pdf"), dpi=300)
