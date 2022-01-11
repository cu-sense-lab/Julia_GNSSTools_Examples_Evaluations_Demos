using GNSSTools
using SatelliteToolbox
using ProgressMeter
using Makie, FileIO, Colors, GeometryBasics, GLMakie
using AbstractPlotting, AbstractPlotting.MakieLayout

GM = 3.986004418e14  # m³s⁻²
eop = get_eop();
user_lla = (40.01, -105.2437, 1655);

###########ONE-WEB############
a = Rₑ+1200*1000;  # meters
plane_num = 18
sat_per_plane = 36;
incl = 87.9;  # degrees
ΔΩ = 360/plane_num/2;
###########IRIDIUM############
a = Rₑ+780*1000;  # meters
plane_num = 6;
sat_per_plane = 11;
incl = 86;
ΔΩ = 360/plane_num/2;
###########STARLINK###########
# a = Rₑ+1150*1000;  # meters
# plane_num = 32;
# sat_per_plane = 50;
# incl = 53;
# ΔΩ = 360/plane_num;
a_starlink = Rₑ+1150*1000;  # meters
#############GPS##############
# a = Rₑ+20000*1000;
# plane_num = 6;
# sat_per_plane = 4;
# incl = 56;
# ΔΩ = 360/plane_num;
##############################
orbital_period = 2π*sqrt(a^3 / GM)  # seconds
orbital_period_starlink = 2π*sqrt(a_starlink^3 / GM)  # seconds
# orb_period_fraction = 2
# orb_period_fraction = (60*60*24) / orbital_period
orb_period_fraction = 1
frames_per_second = 25  # frames/second
total_orbit_time = orb_period_fraction*orbital_period
total_orbit_time_starlink = orb_period_fraction*orbital_period_starlink
total_vid_time = 10  # seconds
# Δt = total_orbit_time/total_vid_time/frames_per_second  # seconds
Δt = total_orbit_time_starlink/total_vid_time/frames_per_second  # seconds
t_range = Array(0:Δt:total_orbit_time);
freq = 9650e6;
constellation = define_constellation(a, plane_num, 
                                     sat_per_plane, incl, 
                                     t_range; obs_lla=user_lla,
                                     print_steps=true,
                                     eop=eop,
                                     ΔΩ=ΔΩ,
                                     show_plot=false);
# Loop through each time step and produce a constellation plot to be combined
# into a video later
obs_ecef = GeodetictoECEF(deg2rad(user_lla[1]), deg2rad(user_lla[2]), user_lla[3])
N = length(constellation.t_range)
sat_num = length(constellation.satellites)
# Earth 
u = Array(range(0, 2π, length=200)) .+ (180*pi/180)
v = Array(range(0, 1π, length=200))
u, v = meshgrid(u, v)
E_x = Rₑ.*cos.(u).*sin.(v)
E_y = Rₑ.*sin.(u).*sin.(v)
E_z = Rₑ.*cos.(v)

# earth = load(download("https://svs.gsfc.nasa.gov/vis/a000000/a002900/a002915/bluemarble-2048.png"));
earth = load("bluemarble-2048.png");

# function make_plot(i, scene=Scene(resolution=(1000, 1000), color=:black);
#                    sat2user_vectors=true)
#     pt_size = 9e8
#     strokewidth = 4
#     linewidth = 3
#     ϕ = deg2rad(270)
#     θ = deg2rad(65)
#     scaler = 3
#     camera_pos = (scaler*constellation.a).*[cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]
#     light = Vec{3,Float32}[[1.0, 1.0, 1.0], [0.1, 0.1, 0.1], [0.9, 0.9, 0.9], [0.0, -20.0, 0.0]]
#     to_eci = rECEFtoECI(ITRF(), TEME(), constellation.epoch+constellation.t_range[i], eop)
#     # Transform Earth and CO from ECEF to ECI frames
#     E_xt = zeros(size(E_x))
#     E_yt = zeros(size(E_y))
#     E_zt = zeros(size(E_z))
#     for j in 1:size(E_x)[1]
#         for k in 1:size(E_x)[2]
#             x, y, z = to_eci*[E_x[j,k], E_y[j,k], E_z[j,k]]
#             E_xt[j,k] = x
#             E_yt[j,k] = y
#             E_zt[j,k] = z
#         end
#     end
#     surface!(E_xt, E_yt, E_zt, color=earth, shading=false, show_axis=false)
#     obs_eci = to_eci*[obs_ecef[1], obs_ecef[2], obs_ecef[3]]
#     for j in 1:sat_num
#         r_ecef = constellation.satellites[j].r_ecef[i,:]
#         r_eci = kepler_to_sv(constellation.satellites[j].orbit[i]).r
#         sat_range, az, el = calcelevation(r_ecef, user_lla)
#         scatter!([r_eci[1]], [r_eci[2]], [r_eci[3]],
#                 strokewidth=strokewidth, strokecolor=:grey)
#         if el > 10
#             lines!([obs_eci[1], r_eci[1]], 
#                    [obs_eci[2], r_eci[2]], 
#                    [obs_eci[3], r_eci[3]], color=:red, linewidth=linewidth)
#         end
#     end
#     for j in 1:sat_num
#         if (j%sat_per_plane) == (sat_per_plane-1)
#             r_ecis = zeros(size(constellation.satellites[j].r_ecef))
#             for k in 1:size(r_ecis)[1]
#                 # r_ecis[k,:] = to_eci*constellation.satellites[j].r_ecef[k,:]
#                 r_ecis[k,:] =  kepler_to_sv(constellation.satellites[j].orbit[k]).r
#             end
#             lines!(r_ecis[:,1], r_ecis[:,2], r_ecis[:,3], color=:grey, linewidth=linewidth)
#         end
#     end
#     update_cam!(scene, camera_pos, 0)
#     return scene
# end

function rotate_earth(i)
    to_eci = rECEFtoECI(ITRF(), TEME(), constellation.epoch+constellation.t_range[i], eop)
    # Transform Earth and CO from ECEF to ECI frames
    E_xt = zeros(size(E_x))
    E_yt = zeros(size(E_y))
    E_zt = zeros(size(E_z))
    for j in 1:size(E_x)[1]
        for k in 1:size(E_x)[2]
            x, y, z = to_eci*[E_x[j,k], E_y[j,k], E_z[j,k]]
            E_xt[j,k] = x
            E_yt[j,k] = y
            E_zt[j,k] = z
        end
    end
    return (E_xt, E_yt, E_zt)
end

function get_obj_positions(i)
    to_eci = rECEFtoECI(ITRF(), TEME(), constellation.epoch+constellation.t_range[i], eop)
    obs_eci = to_eci*obs_ecef
    sat_pos = zeros(length(constellation.satellites), 3)
    planes = []
    visible = []
    visible_sat = []
    for j in 1:length(constellation.satellites)
        r_ecef = constellation.satellites[j].r_ecef[i,:]
        r_eci = kepler_to_sv(constellation.satellites[j].orbit[i]).r
        sat_pos[j,:] = r_eci
        sat_range, az, el = calcelevation(r_ecef, user_lla)
        if el > 10
            push!(visible, [[obs_eci[1], r_eci[1]], [obs_eci[2], r_eci[2]], [obs_eci[3], r_eci[3]]])
            push!(visible_sat, r_eci)
        end
        if (constellation.satellites[j].id%sat_per_plane) == 1
            r_ecis = zeros(size(constellation.satellites[j].r_ecef))
            for k in 1:size(r_ecis)[1]
                r_ecis[k,:] =  kepler_to_sv(constellation.satellites[j].orbit[k]).r
            end
            push!(planes, r_ecis)
        end
    end
    visible_sats = zeros(length(sat_pos),3)
    for j in 1:size(sat_pos)[1]
        if j > length(visible_sat)
            visible_sats[j,:] = [NaN, NaN, NaN]
        else
            visible_sats[j,:] = visible_sat[j]
        end
    end
    return (sat_pos, planes, visible, visible_sats)
end

function initialize_figure(i; show_vectors=true, visible_color=missing)
    ϕ = deg2rad(180)
    θ = deg2rad(70)
    scaler = 3
    camera_pos = (scaler*constellation.a).*[cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]
    light_ϕ = deg2rad(270)
    light_θ = deg2rad(90)
    light_scaler = 149597870700/Rₑ
    light_camera_pos = (light_scaler*constellation.a).*[cos(light_ϕ)*sin(light_θ), sin(light_ϕ)*sin(light_θ), cos(light_θ)]
    lightpos = Vec3f0(light_camera_pos...)
    pt_size = 9e8
    strokewidth = 3
    linewidth = 2
    linestyle = nothing
    pt_color = :lightgrey
    line_color = :lightgrey
    vector_color = :lightblue
    background_color = :white
    if ismissing(visible_color)
        visible_color = pt_color
    end
    E_xt, E_yt, E_zt = rotate_earth(i)
    sat_pos, planes, visible, visible_sats = get_obj_positions(i)
    scene = Scene(resolution=(1000, 1000), backgroundcolor=background_color)
    surf = surface!(scene, E_xt, E_yt, E_zt, shading=true, color=earth, 
                    show_axis=false, lightposition=lightpos, 
                    colormap=(:white, :white),
                    ambient = Vec3f0(0.5, 0.5, 0.5), 
                    diffuse = Vec3f0(0.4, 0.4, 0.2),
                    specular = Vec3f0(0.5, 0.5, 0.5), 
                    shininess = 0.5f0)
    scatter_sat = scatter!(scene, sat_pos[:,1], sat_pos[:,2], 
    sat_pos[:,3],
                           strokewidth=strokewidth, strokecolor=pt_color)
    scatter_visible_sat = scatter!(scene, 
                                   visible_sats[:,1], 
                                   visible_sats[:,2],
                                   visible_sats[:,3], 
                                   strokewidth=strokewidth, 
                                   strokecolor=visible_color)
    for j in 1:length(planes)
        line_planes = lines!(scene, planes[j][:,1], planes[j][:,2],  
                             planes[j][:,3],
                             color=line_color, linewidth=linewidth,
                             linestyle=linestyle)
    end
    if show_vectors
        val = [NaN, NaN]
        for j in 1:ceil(Int, length(constellation.satellites)/2)
            line_visible = lines!(scene, val, val,  val,
                                  color=vector_color, linewidth=linewidth)
        end
        for j in 1:length(visible)
            scene.plots[j+length(planes)+3][1] = visible[j][1]
            scene.plots[j+length(planes)+3][2] = visible[j][2]
            scene.plots[j+length(planes)+3][3] = visible[j][3]
        end
    end
    update_cam!(scene, camera_pos, 0)
    return scene
end

function modify_figure!(scene, i; show_vectors=true)
    ϕ = deg2rad(180)
    θ = deg2rad(70)
    scaler = 3
    camera_pos = (scaler*constellation.a).*[cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]
    E_xt, E_yt, E_zt = rotate_earth(i)
    sat_pos, planes, visible, visible_sats = get_obj_positions(i)
    scene.plots[1][1] = E_xt
    scene.plots[1][2] = E_yt
    scene.plots[1][3] = E_zt
    scene.plots[2][1] = sat_pos[:,1]
    scene.plots[2][2] = sat_pos[:,2]
    scene.plots[2][3] = sat_pos[:,3]
    scene.plots[3][1] = visible_sats[:,1]
    scene.plots[3][2] = visible_sats[:,2]
    scene.plots[3][3] = visible_sats[:,3]
    for j in 1:length(planes)
        scene.plots[j+3][1] = planes[j][:,1]
        scene.plots[j+3][2] = planes[j][:,2]
        scene.plots[j+3][3] = planes[j][:,3]
    end
    if show_vectors
        val = [NaN, NaN]
        for j in 1:(length(scene.plots)-3-length(planes))
            scene.plots[j+length(planes)+3][1] = val
            scene.plots[j+length(planes)+3][2] = val
            scene.plots[j+length(planes)+3][3] = val

        end
        for j in 1:length(visible)  
            scene.plots[j+length(planes)+3][1] = visible[j][1]
            scene.plots[j+length(planes)+3][2] = visible[j][2]
            scene.plots[j+length(planes)+3][3] = visible[j][3]
        end
    end
    update_cam!(scene, camera_pos, 0)
    return scene
end

function make_plots(N; show_vectors=true, visible_color=missing)
    scene = initialize_figure(1; show_vectors=show_vectors, 
                              visible_color=visible_color)
    p = Progress(N, 1, "Generating $N frames...")
    compression = 10
    record(scene, "gps_constellation_black_compressed.mp4", 1:N; sleep=false, 
           compression=compression, framerate=frames_per_second) do i
        modify_figure!(scene, i; show_vectors=show_vectors)
        next!(p)
    end
end

function zoom_earth(N=11*frames_per_second, start_i=9*frames_per_second, 
                    end_i=10*frames_per_second; 
                    show_vectors=true, visible_color=missing)
    a1 = Rₑ+20000*1000
    a2 = Rₑ+1150*1000;
    ϕ = deg2rad(180)
    θ = deg2rad(70)
    scaler = 3
    camera_pos = (scaler*a1).*[cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]
    light_ϕ = deg2rad(270)
    light_θ = deg2rad(90)
    light_scaler = 149597870700/Rₑ
    light_camera_pos = (light_scaler*a1).*[cos(light_ϕ)*sin(light_θ), sin(light_ϕ)*sin(light_θ), cos(light_θ)]
    lightpos = Vec3f0(light_camera_pos...)
    E_xt, E_yt, E_zt = rotate_earth(1)
    scene = Scene(resolution=(1000, 1000), backgroundcolor=:black)
    surf = surface!(scene, E_xt, E_yt, E_zt, shading=true, color=earth, 
                    show_axis=false, lightposition=lightpos, 
                    colormap=(:white, :white),
                    ambient = Vec3f0(0.5, 0.5, 0.5), 
                    diffuse = Vec3f0(0.4, 0.4, 0.2),
                    specular = Vec3f0(0.5, 0.5, 0.5), 
                    shininess = 0.5f0)
    as = fill(a1, N)
    Δas = (a2 - a1)
    for i in start_i:end_i
        as[i] = a1 - Δas*(cos((i - start_i) * π / (end_i - start_i))/2 - 0.5)
    end
    as[end_i+1:end] .= a2
    p = Progress(N, 1, "Generating $N frames...")
    compression = 0
    record(scene, "zoom_earth.mp4", 1:N; sleep=false, 
           compression=compression, framerate=frames_per_second) do i
        E_xt, E_yt, E_zt = rotate_earth(i)
        scene.plots[1][1] = E_xt
        scene.plots[1][2] = E_yt
        scene.plots[1][3] = E_zt
        camera_pos = (scaler*as[i]).*[cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]
        update_cam!(scene, camera_pos, 0)
        next!(p)
    end
end

function make_hist_plots()
    # GPS
    a = Rₑ+20000*1000;
    plane_num = 6;
    sat_per_plane = 4;
    incl = 56;
    ΔΩ = 360/plane_num;
    orbital_period = 2π*sqrt(a^3 / GM)  # seconds
    # orb_period_fraction = 2
    # orb_period_fraction = (60*60*24) / orbital_period
    orb_period_fraction = 10
    frames_per_second = 1  # frames/second
    total_orbit_time = orb_period_fraction*orbital_period
    total_vid_time = total_orbit_time
    # Δt = total_orbit_time/total_vid_time/frames_per_second  # seconds
    Δt = total_orbit_time/total_vid_time/frames_per_second  # seconds
    t_range = Array(0:Δt:total_orbit_time);
    # freq = 9650e6;
    freq = L1_freq
    doppler, doppler_rate, doppler_bounds, 
    doppler_rate_bounds = doppler_distribution(a, plane_num, 
                                               sat_per_plane, incl, 
                                               t_range, user_lla, 
                                               freq, show_plot=false,
                                               print_steps=true, 
                                               eop=eop,
                                               ΔΩ=ΔΩ);
    # Make plot
    scene, layout = layoutscene(resolution=(840, 360))
    ax1 = layout[1, 1] = LAxis(scene, xlabel="Doppler Frequency (kHz)")
    hideydecorations!(ax1, label=true, ticklabels=true, ticks=true, grid=true)
    hidexdecorations!(ax1, label=false, ticklabels=false, ticks=false, grid=true)
    barplot!(ax1, doppler.edges[1][2:end]./1000, doppler.weights, color=:black, 
            strokewidth=1, strokecolor=:black, show_axis=true)
    save("gps_hist_doppler.png", scene)
    scene, layout = layoutscene(resolution=(840, 360))
    ax1 = layout[1, 1] = LAxis(scene, xlabel = "Doppler Frequency Rate (Hz/s)")    
    hideydecorations!(ax1, label=true, ticklabels=true, ticks=true, grid=true)
    hidexdecorations!(ax1, label=false, ticklabels=false, ticks=false, grid=true)
    barplot!(ax1, doppler_rate.edges[1][2:end], doppler_rate.weights, color=:black, 
            strokewidth=1, strokecolor=:black, show_axis=true)
    save("gps_hist_doppler_rate.png", scene)
    # Starlink
    a = Rₑ+1150*1000;  # meters
    plane_num = 32;
    sat_per_plane = 50;
    incl = 53;
    ΔΩ = 360/plane_num;
    orbital_period = 2π*sqrt(a^3 / GM)  # seconds
    # orb_period_fraction = 2
    # orb_period_fraction = (60*60*24) / orbital_period
    orb_period_fraction = 10
    frames_per_second = 1  # frames/second
    total_orbit_time = orb_period_fraction*orbital_period
    total_vid_time = total_orbit_time
    # Δt = total_orbit_time/total_vid_time/frames_per_second  # seconds
    Δt = total_orbit_time/total_vid_time/frames_per_second  # seconds
    t_range = Array(0:Δt:total_orbit_time);
    # freq = 9650e6;
    freq = L1_freq
    doppler, doppler_rate, doppler_bounds, 
    doppler_rate_bounds = doppler_distribution(a, plane_num, 
                                               sat_per_plane, incl, 
                                               t_range, user_lla, 
                                               freq, show_plot=false,
                                               print_steps=true, 
                                               eop=eop,
                                               ΔΩ=ΔΩ);
    # Make plot
    scene, layout = layoutscene(resolution=(840, 360))
    ax1 = layout[1, 1] = LAxis(scene, xlabel="Doppler Frequency (kHz)")
    hideydecorations!(ax1, label=true, ticklabels=true, ticks=true, grid=true)
    hidexdecorations!(ax1, label=false, ticklabels=false, ticks=false, grid=true)
    barplot!(ax1, doppler.edges[1][2:end]./1000, doppler.weights, color=:black, 
            strokewidth=1, strokecolor=:black, show_axis=true)
    save("starlink_hist_doppler.png", scene)
    scene, layout = layoutscene(resolution=(840, 360))
    ax1 = layout[1, 1] = LAxis(scene, xlabel = "Doppler Frequency Rate (Hz/s)")    
    hideydecorations!(ax1, label=true, ticklabels=true, ticks=true, grid=true)
    hidexdecorations!(ax1, label=false, ticklabels=false, ticks=false, grid=true)
    barplot!(ax1, doppler_rate.edges[1][2:end], doppler_rate.weights, color=:black, 
            strokewidth=1, strokecolor=:black, show_axis=true)
    ax1.xticks = -200:50:0
    save("starlink_hist_doppler_rate.png", scene)
end
# run(`ffmpeg -r $frames_per_second -f image2 -i constellation_figures/constellation_%0$(pad_num)d.png -y constellation.mp4`)
run(`ffmpeg -i LabVideo_Sergei.mp4 -filter_complex "[0:v] palettegen" -y palette.png`)
run(`ffmpeg -i LabVideo_Sergei.mp4 -i palette.png -filter_complex "[0:v][1:v] paletteuse" -y LabVideo_Sergei.gif`)


# L1 C/A signal undergoing extreme doppler
using PyPlot
pygui(false)
t_length = 2
sig_freq = 9650e6;
sigtype = define_l1ca_code_type(t_length; sig_freq=sig_freq);
doppler_t = Array(0:0.001:t_length+0.1);
f_d = 150*1e3
fd_rate = -1000
fdd_rate = 175
fddd_rate = 0
doppler_curve = f_d .+ fd_rate .* doppler_t .+ fdd_rate .* doppler_t.^2 .+ fddd_rate .* doppler_t.^3;
signal = definesignal(sigtype, 5e6, t_length);
generatesignal!(signal; doppler_curve=doppler_curve, doppler_t=doppler_t);
acqresults, trackresults = process(signal, sigtype, 1; fd_center=f_d, show_plot=false);
# Iridium constellation @ X-band
t_range = Array(0:1:12*60*60);
a = Rₑ+780*1000;
plane_num = 6;
sat_per_plane = 11;
incl = 86;
freq = 9650e6;
ΔΩ = 30;
doppler, doppler_rate, doppler_bounds, 
doppler_rate_bounds = doppler_distribution(a, plane_num, 
                                           sat_per_plane, incl, 
                                           t_range, user_lla, 
                                           freq, 
                                           print_steps=true, 
                                           eop=eop,
                                           ΔΩ=ΔΩ;
                                           show_plot=false);
# Make animated plots
figsize = (17,8)
frames_per_second = 25
total_length = 10
N = floor(Int, frames_per_second*t_length)
pad_num = length(string(N))
ΔN = floor(Int, length(trackresults.t)/N)
p = Progress(N, 1, "Generating $N frames...")
for i in 0:N
    idx = i*ΔN
    fig = figure(figsize=figsize)
    ax1 = subplot(221)
    ax1.plot(trackresults.t[1:idx], trackresults.n0[1:idx].%1023, "k")
    xlim([trackresults.t[1], trackresults.t[end]])
    ylim([minimum(trackresults.n0.%1023), maximum(trackresults.n0.%1023)])
    title("Code Phase", fontsize=20)
    xticks([])
    yticks([])
    ax2 = subplot(222)
    ax2.plot(trackresults.t[1:idx], trackresults.dphi_meas[1:idx], "k")
    xlim([trackresults.t[1], trackresults.t[end]])
    ylim([minimum(trackresults.dphi_meas), maximum(trackresults.dphi_meas)])
    title("Phase Residuals", fontsize=20)
    xticks([])
    yticks([])
    ax3 = subplot(223)
    ax3.plot(trackresults.t[1:idx], trackresults.fds[1:idx], "k")
    xlim([trackresults.t[1], trackresults.t[end]])
    ylim([minimum(trackresults.fds), maximum(trackresults.fds)])
    title("Doppler", fontsize=20)
    xticks([])
    yticks([])
    ax4 = subplot(224)
    ax4.plot(trackresults.t[1:idx], real.(trackresults.ZP[1:idx]), "k")
    xlim([trackresults.t[1], trackresults.t[end]])
    ylim([minimum(real.(trackresults.ZP)), maximum(real.(trackresults.ZP))])
    # ax4.plot(trackresults.t, trackresults.data_bits, "k", linestyle="-")
    title("Navigation Data", fontsize=20)
    xticks([])
    yticks([])
    subplots_adjust(top=0.96, bottom=0.004, right=0.998, left=0.001, 
                    hspace=0.2, wspace=0.1)
    # margins(0,0,0)
    # gca().xaxis.set_major_locator(PyPlot.plt.NullLocator())
    # gca().yaxis.set_major_locator(PyPlot.plt.NullLocator())
    if i == N
        num = string(0, pad=pad_num)
    else
        num = string(i+1, pad=pad_num)
    end
    savefig("data_figures/data_$num.png")
    close(fig)
    next!(p)
end
run(`ffmpeg -r $frames_per_second -f image2 -i data_figures/data_%0$(pad_num)d.png -y data.mp4`)
run(`ffmpeg -i data.mp4 -filter_complex "[0:v] palettegen" -y palette1.png`)
run(`ffmpeg -i data.mp4 -i palette1.png -filter_complex "[0:v][1:v] paletteuse" -y data.gif`)