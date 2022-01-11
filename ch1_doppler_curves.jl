using GNSSTools
using JLD
using PyPlot
pygui(true)

GM = 3.986004418e14  # m³s⁻²
eop = get_eop();

# N_rand = 25
# lats = rand(-90:0.001:90, N_rand)  # degrees
# longs = rand(-180:0.001:180, N_rand)  # degrees
# heights = fill(100, N_rand)  # meters
# llas = (lats, longs, heights)
# user_lla = llas
freq = L1_freq
# user_lla = (40.01, -105.2437, 1655);

file_name = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/ch1_doppler_curves_user_lla.jld"
# file_name = "ch1_doppler_curves_user_lla.jld"
user_lla = load(file_name, "user_lla")


Δt = 1  # seconds
total_time = 1*24*60*60  # seconds
t_range = Array(0:Δt:total_time);
min_elevation = 10  # degrees

# GPS
a = 26_560e3;  # m
plane_num = 6;
sat_per_plane = 4;
incl = 56;
ΔΩ = 360/plane_num;
Δf_per_plane = 360/sat_per_plane/2;  # degrees
orbital_period = 2π*sqrt(a^3 / GM)  # seconds
doppler_gps, doppler_rate_gps, doppler_bounds_gps, 
doppler_rate_bounds_gps, raw_dopplers_gps = doppler_distribution(a, plane_num, 
                                           sat_per_plane, incl, 
                                           t_range, user_lla, 
                                           freq, show_plot=false,
                                           print_steps=true, 
                                           eop=eop,
                                           min_elevation=min_elevation,
                                           ΔΩ=ΔΩ,
                                           Δf_per_plane=Δf_per_plane,
										   return_constellation=false,
                                           return_raw_dopplers=true);
doppler_raw_gps, doppler_rate_raw_gps, 
elevations_gps, max_elevations_gps = raw_dopplers_gps

# Iridium
a = Rₑ+780*1000;  # meters
plane_num = 6;
sat_per_plane = 11;
incl = 86;
ΔΩ = 360/plane_num/2;
Δf_per_plane = 360/sat_per_plane/2;  # degrees
orbital_period = 2π*sqrt(a^3 / GM)  # seconds
doppler_iridium, doppler_rate_iridium, doppler_bounds_iridium, 
doppler_rate_bounds_iridium, raw_dopplers_iridium = doppler_distribution(a, plane_num, 
                                            sat_per_plane, incl, 
                                            t_range, user_lla, 
                                            freq, show_plot=false,
                                            print_steps=true, 
                                            eop=eop,
                                            ΔΩ=ΔΩ,
                                            Δf_per_plane=Δf_per_plane,
                                            return_raw_dopplers=true);
doppler_raw_iridium, doppler_rate_raw_iridium, 
elevations_iridium, max_elevations_iridium = raw_dopplers_iridium

# Starlink
a = Rₑ+550*1000;  # meters
plane_num = 24;
sat_per_plane = 66;
incl = 53;
ΔΩ = 360/plane_num;
Δf_per_plane = 360/sat_per_plane/2;  # degrees
orbital_period = 2π*sqrt(a^3 / GM)  # seconds
doppler_starlink, doppler_rate_starlink, doppler_bounds_starlink, 
doppler_rate_bounds_starlink, raw_dopplers_starlink = doppler_distribution(a, plane_num, 
                                            sat_per_plane, incl, 
                                            t_range, user_lla, 
                                            freq, show_plot=false,
                                            print_steps=true, 
                                            eop=eop,
                                            ΔΩ=ΔΩ,
                                            Δf_per_plane=Δf_per_plane,
                                            return_raw_dopplers=true);
doppler_raw_starlink, doppler_rate_raw_starlink, 
elevations_starlink, max_elevations_starlink = raw_dopplers_starlink

# OneWeb
a = Rₑ+1200*1000;  # meters
plane_num = 18;
sat_per_plane = 36;
incl = 87.9;  # degrees
ΔΩ = 360/plane_num/2;
Δf_per_plane = 360/sat_per_plane/2;  # degrees
orbital_period = 2π*sqrt(a^3 / GM)  # seconds
doppler_oneweb, doppler_rate_oneweb, doppler_bounds_oneweb, 
doppler_rate_bounds_oneweb, raw_dopplers_oneweb = doppler_distribution(a, plane_num, 
                                            sat_per_plane, incl, 
                                            t_range, user_lla, 
                                            freq, show_plot=false,
                                            print_steps=true, 
                                            eop=eop,
                                            ΔΩ=ΔΩ,
                                            Δf_per_plane=Δf_per_plane,
                                            return_raw_dopplers=true);
doppler_raw_oneweb, doppler_rate_raw_oneweb, 
elevations_oneweb, max_elevations_oneweb = raw_dopplers_oneweb

save("ch1_doppler_curves.jld", 
     "doppler_gps", doppler_gps,
	 "doppler_iridium", doppler_iridium,
	 "doppler_starlink", doppler_starlink,
	 "doppler_oneweb", doppler_oneweb,
	 "doppler_rate_gps", doppler_rate_gps,
	 "doppler_rate_iridium", doppler_rate_iridium,
	 "doppler_rate_starlink", doppler_rate_starlink,
	 "doppler_rate_oneweb", doppler_rate_oneweb,
	 "doppler_bounds_gps", doppler_bounds_gps,
	 "doppler_bounds_iridium", doppler_bounds_iridium,
	 "doppler_bounds_starlink", doppler_bounds_starlink,
	 "doppler_bounds_oneweb", doppler_bounds_oneweb,
	 "doppler_raw_gps", doppler_raw_gps,
	 "doppler_raw_iridium", doppler_raw_iridium,
	 "doppler_raw_starlink", doppler_raw_starlink,
	 "doppler_raw_oneweb", doppler_raw_oneweb,
	 "doppler_rate_raw_gps", doppler_rate_raw_gps,
	 "doppler_rate_raw_iridium", doppler_rate_raw_iridium,
	 "doppler_rate_raw_starlink", doppler_rate_raw_starlink,
	 "doppler_rate_raw_oneweb", doppler_rate_raw_oneweb,
	 "freq", freq,
	 "user_lla", user_lla,
	 "t_range", t_range,
	 "min_elevation", min_elevation,
	 "elevations_gps", elevations_gps,
	 "elevations_iridium", elevations_iridium,
	 "elevations_starlink", elevations_starlink,
	 "elevations_oneweb", elevations_oneweb,
	 "max_elevations_gps", max_elevations_gps,
	 "max_elevations_iridium", max_elevations_iridium,
	 "max_elevations_starlink", max_elevations_starlink,
	 "max_elevations_oneweb", max_elevations_oneweb)


# file_name = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/ch1_doppler_curves.jld"
# data = load(file_name, 
#             "doppler_gps",
# 	       "doppler_iridium",
# 	       "doppler_starlink",
# 	       "doppler_oneweb",
# 	       "doppler_rate_gps",
# 	       "doppler_rate_iridium",
# 	       "doppler_rate_starlink",
# 	       "doppler_rate_oneweb",
# 	       "doppler_bounds_gps",
# 	       "doppler_bounds_iridium",
# 	       "doppler_bounds_starlink",
# 	       "doppler_bounds_oneweb",
# 	       "doppler_raw_gps",
# 	       "doppler_raw_iridium",
# 	       "doppler_raw_starlink",
# 	       "doppler_raw_oneweb",
# 	       "doppler_rate_raw_gps",
# 	       "doppler_rate_raw_iridium",
# 	       "doppler_rate_raw_starlink",
# 	       "doppler_rate_raw_oneweb",
# 	       "freq",
# 	       "user_lla",
# 	       "t_range",
# 	       "min_elevation")


# doppler_gps, doppler_iridium, doppler_starlink, doppler_oneweb, 
# doppler_rate_gps, doppler_rate_iridium, doppler_rate_starlink, doppler_rate_oneweb,
# doppler_bounds_gps, doppler_bounds_iridium, doppler_bounds_starlink, doppler_bounds_oneweb,
# doppler_raw_gps, doppler_raw_iridium, doppler_raw_starlink, doppler_raw_oneweb,
# doppler_rate_raw_gps, doppler_rate_raw_iridium, doppler_rate_raw_starlink, doppler_rate_raw_oneweb,
# freq, user_lla, t_range, min_elevation = data 

# # Make plots
# fig = figure(figsize=(8.5,9))
# # GPS
# ax1 = fig.add_subplot(4,2,1)
# ax1.hist(Array(doppler_gps.edges[1])[2:end]./1000,
#          Array(doppler_gps.edges[1])./1000, weights=doppler_gps.weights,
#          color="k");
# xlabel("Doppler (kHz)")
# ylabel("Probability")
# title("GPS")
# ax2 = fig.add_subplot(4,2,2)
# ax2.hist(Array(doppler_rate_gps.edges[1])[2:end],
#          Array(doppler_rate_gps.edges[1]), weights=doppler_rate_gps.weights,
#          color="k");
# xlabel("Doppler Rate (Hz/s)")
# ylabel("Probability")
# title("GPS")
# # Iridium
# ax3 = fig.add_subplot(4,2,3)
# ax3.hist(Array(doppler_iridium.edges[1])[2:end]./1000,
#          Array(doppler_iridium.edges[1])./1000, weights=doppler_iridium.weights,
#          color="k");
# xlabel("Doppler (kHz)")
# ylabel("Probability")
# title("Iridium")
# ax4 = fig.add_subplot(4,2,4)
# ax4.hist(Array(doppler_rate_iridium.edges[1])[2:end],
#          Array(doppler_rate_iridium.edges[1]), weights=doppler_rate_iridium.weights,
#          color="k");
# xlabel("Doppler Rate (Hz/s)")
# ylabel("Probability")
# title("Iridium")
# # Starlink
# ax5 = fig.add_subplot(4,2,5)
# ax5.hist(Array(doppler_starlink.edges[1])[2:end]./1000,
#          Array(doppler_starlink.edges[1])./1000, weights=doppler_starlink.weights,
#          color="k");
# xlabel("Doppler (kHz)")
# ylabel("Probability")
# title("Starlink")
# ax6 = fig.add_subplot(4,2,6)
# ax6.hist(Array(doppler_rate_starlink.edges[1])[2:end],
#          Array(doppler_rate_starlink.edges[1]), weights=doppler_rate_starlink.weights,
#          color="k");
# xlabel("Doppler Rate (Hz/s)")
# ylabel("Probability")
# title("Starlink")
# # OneWeb
# ax7 = fig.add_subplot(4,2,7)
# ax7.hist(Array(doppler_oneweb.edges[1])[2:end]./1000,
#          Array(doppler_oneweb.edges[1])./1000, weights=doppler_oneweb.weights,
#          color="k");
# xlabel("Doppler (kHz)")
# ylabel("Probability")
# title("OneWeb")
# ax8 = fig.add_subplot(4,2,8)
# ax8.hist(Array(doppler_rate_oneweb.edges[1])[2:end],
#          Array(doppler_rate_oneweb.edges[1]), weights=doppler_rate_oneweb.weights,
#          color="k");
# xlabel("Doppler Rate (Hz/s)")
# ylabel("Probability")
# title("OneWeb")

# subplots_adjust(hspace=0.7, wspace=0.35, top=0.93, left=0.1, right=0.9)
# savefig("/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/figures/ch1_doppler_distributions.svg", dpi=300)


# using PyCall
# basemap = pyimport("mpl_toolkits.basemap")
# # lon_0 is central longitude of projection.
# # resolution = 'c' means use crude resolution coastlines.
# fig = figure(figsize=(7.5, 2.5))
# ax1 = fig.add_subplot(1,1,1)
# m = basemap.Basemap(projection="robin", lon_0=0, resolution="c")
# m.drawcoastlines(color="black", linewidth=0.25)
# m.fillcontinents(color="#69b2a2",lake_color="#A6CAE0")
# # draw parallels and meridians.
# m.drawparallels(range(-90., 120., step=30.), linewidth=0.5, color="dimgrey")
# m.drawmeridians(range(0., 360., step=60.), linewidth=0.5, color="dimgrey")
# m.drawmapboundary(fill_color="#A6CAE0")
# # m.fillcontinents(color="grey", alpha=0.3)
# color = "tomato"
# marker = "."
# markersize = 5
# N = min(length(user_lla[1]), length(user_lla[2]), length(user_lla[3]))
# for i in 1:N
# 	lat = user_lla[1][i]
# 	lon	= user_lla[2][i]
# 	xpt, ypt = m(lon,lat)
# 	m.plot(xpt, ypt, color=color, marker=marker, markersize=markersize)
# end
# # title("$N Random Observation Locations")
# subplots_adjust(bottom=0.01, left=0.1, right=0.9, top=0.99)
# savefig("/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/figures/ch1_doppler_curves_random_location_map.pdf",
#         dpi=300)
