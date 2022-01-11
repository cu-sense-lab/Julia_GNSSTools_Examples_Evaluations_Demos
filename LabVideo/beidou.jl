using GNSSTools
using SatelliteToolbox
using ProgressMeter
using Makie
using GLMakie
using FileIO
# using GLMakie: Scene, surface!, scatter!, lines!, update_cam!, save
using Colors, GeometryBasics
# using AbstractPlotting, AbstractPlotting.MakieLayout


"""
Constellation info from 
[www.beidou.gov.cn](http://en.beidou.gov.cn/SYSTEMS/Officialdocument/202110/P020211014595952404052.pdf)
"""


GM = 3.986004418e14  # m³s⁻²
eop = get_eop();
scale_a = Rₑ + 35786*1000

Ω₀ = -30
a = Rₑ+21528*1000;  # meters
plane_num = 3;
sat_per_plane = 8;
incl = 55;  # degrees
ΔΩ = 360/plane_num;  # degrees
Δf_per_plane = 0;  # degrees
orbital_period = 2π*sqrt(a^3 / GM)  # seconds


Δt = 1 # seconds
t_range = Array(0:Δt:orbital_period);  # seconds


constellation = define_constellation(a, plane_num, 
                                     sat_per_plane, incl, 
                                     t_range; obs_lla=missing,
                                     print_steps=true,
                                     eop=eop,
                                     ΔΩ=ΔΩ,
                                     show_plot=false,
                                     Δf_per_plane,
                                     Ω₀=Ω₀);

# Earth sphere
u = Array(range(0, 2π, length=200)) .- (55*pi/180);
v = Array(range(0, 1π, length=200));
u, v = meshgrid(u, v);
E_x = Rₑ.*cos.(u).*sin.(v);
E_y = Rₑ.*sin.(u).*sin.(v);
E_z = Rₑ.*cos.(v);


# Earth imagery
# earth = load(download("https://svs.gsfc.nasa.gov/vis/a000000/a002900/a002915/bluemarble-2048.png"));
earth = load("bluemarble-2048.png");


# View angle
ϕ = deg2rad(120)
θ = deg2rad(80)
scaler = 3
camera_pos = (scaler*scale_a).*[cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]


# Parameters
light_ϕ = deg2rad(270)
light_θ = deg2rad(90)
light_scaler = 149597870700/Rₑ
light_camera_pos = (light_scaler*scale_a).*[cos(light_ϕ)*sin(light_θ), sin(light_ϕ)*sin(light_θ), cos(light_θ)]
lightpos = Vec3f0(light_camera_pos...)
resolution = (1000, 1000)
save_resolution_scale = 2
pt_scale = 2
pt_size = (save_resolution_scale*pt_scale) * 9e8
strokewidth = (save_resolution_scale*pt_scale) * 3
linewidth = (save_resolution_scale*pt_scale) * 2
linestyle = nothing
pt_color = :lightgrey
line_color = :lightgrey
background_color = :white


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


function get_obj_positions(constellation)
    i = 1
    sat_num = length(constellation.satellites)
    to_eci = rECEFtoECI(ITRF(), TEME(), constellation.epoch+constellation.t_range[i], eop)
    sat_pos = zeros(length(constellation.satellites), 3)
    planes = []
    for j in 1:length(constellation.satellites)
        t = constellation.satellites[j].t[i]
        to_eci = rECEFtoECI(ITRF(), TEME(), t, eop)
        r_ecef = constellation.satellites[j].r_ecef[i,:]
        # r_eci = kepler_to_sv(constellation.satellites[j].orbit[i]).r
        r_eci = to_eci*r_ecef
        sat_pos[j,:] = r_eci
        if (constellation.satellites[j].id%sat_per_plane) == 1
            r_ecis = zeros(size(constellation.satellites[j].r_ecef))
            for k in 1:size(r_ecis)[1]
                t = constellation.satellites[j].t[k]
                to_eci = rECEFtoECI(ITRF(), TEME(), t, eop)
                # r_ecis[k,:] =  kepler_to_sv(constellation.satellites[j].orbit[k]).r
                r_ecis[k,:] = to_eci*constellation.satellites[j].r_ecef[k,:]
            end
            push!(planes, r_ecis)
        end
    end
    return (sat_pos, planes)
end

# Create scene
scene = Scene(resolution=resolution, backgroundcolor=background_color)

sat_pos, planes = get_obj_positions(constellation)
E_xt, E_yt, E_zt = rotate_earth(1)
surf = surface!(scene, E_xt, E_yt, E_zt, shading=true, color=earth, 
                show_axis=false, #lightposition=lightpos, 
                colormap=(:white, :white),
                ambient = Vec3f0(0.5, 0.5, 0.5), 
                diffuse = Vec3f0(0.4, 0.4, 0.2),
                specular = Vec3f0(0.5, 0.5, 0.5), 
                shininess = 0.5f0)
scatter_sat = scatter!(scene, sat_pos[:,1], sat_pos[:,2], sat_pos[:,3],
                       strokewidth=strokewidth, strokecolor=pt_color)
for j in 1:length(planes)
    line_planes = lines!(scene, planes[j][:,1], planes[j][:,2],  
                         planes[j][:,3],
                         color=line_color, linewidth=linewidth,
                         linestyle=linestyle)
end


# GEO satellites
# a = 42164e3
a = Rₑ + 35786*1000
# Ωs = [158.75, 180.0, 210.5, 240.0, 260.0]  # degrees
Ωs = [80.0, 110.5, 140.0]  # degrees
eccentricity = 0.0
inclination = 0.0
true_anomaly = 0.0
arg_of_perigee = 0.0
orbital_period_geo = 2π*sqrt(a^3 / GM)  # seconds
t_range_geo = Array(0:Δt:orbital_period_geo);  # seconds
geo_pt_color = :lightblue
geo_line_color = geo_pt_color
scatter_geo = zeros(length(Ωs), 3)
for i in 1:length(Ωs)
    Ω = deg2rad(Ωs[i])
    init_orbit = init_orbit_propagator(Val(:twobody), 0., a, eccentricity, 
                                       inclination, Ω, arg_of_perigee, 
                                       true_anomaly)
    r, v = propagate_to_epoch!(init_orbit, t_range_geo)
    scatter_geo[i,:] = r[1]  # ECI TEME
    if i == 1
        plane_geo = Array{Float64}(undef, length(r), 3)
        for j in 1:length(r)
            plane_geo[j,:] = r[j]
        end
        line_planes_geo = lines!(scene, plane_geo[:,1], plane_geo[:,2],  
                                 plane_geo[:,3],
                                 color=geo_line_color, linewidth=linewidth,
                                 linestyle=linestyle)
    end
end
scatter_sat = scatter!(scene, scatter_geo[:,1], scatter_geo[:,2], scatter_geo[:,3],
                       strokewidth=strokewidth, strokecolor=geo_pt_color)


# IGSO satellites
# a = 42164e3
a = Rₑ + 35786*1000
# Ωs = [218.0, 98.0, 338.0]  # degrees
Ωs = [0, 1, 2] .* ΔΩ .+ Ω₀
eccentricity = 0.0
inclination = deg2rad(55.0)
# true_anomaly = [218.0, 98.0, 338.0]
true_anomaly = deepcopy(Ωs)
arg_of_perigee = 0.0
orbital_period_igso = 2π*sqrt(a^3 / GM)  # seconds
t_range_igso = Array(0:Δt:orbital_period_igso);  # seconds
igso_pt_color = :lightpink
igso_line_color = igso_pt_color
scatter_igso = zeros(length(Ωs), 3)
for i in 1:length(Ωs)
    Ω = deg2rad(Ωs[i])
    f = deg2rad(true_anomaly[i])
    init_orbit = init_orbit_propagator(Val(:twobody), 0., a, eccentricity, 
                                       inclination, Ω, arg_of_perigee, 
                                       f)
    r, v = propagate_to_epoch!(init_orbit, t_range_geo)
    scatter_igso[i,:] = r[1]  # ECI TEME
    plane_igso = Array{Float64}(undef, length(r), 3)
    for j in 1:length(r)
        plane_igso[j,:] = r[j]
    end
    line_planes_igso = lines!(scene, plane_igso[:,1], plane_igso[:,2],  
                              plane_igso[:,3],
                              color=igso_line_color, linewidth=linewidth,
                              linestyle=linestyle)
end
scatter_sat = scatter!(scene, scatter_igso[:,1], scatter_igso[:,2], scatter_igso[:,3],
                       strokewidth=strokewidth, strokecolor=igso_pt_color)

display(scene)
update_cam!(scene, cameracontrols(scene), Vec3f0(camera_pos), Vec3f0(0), Vec3f0(0, 0, 1))
scene.center=false; save("beidou.png", scene, resolution=(save_resolution_scale*resolution[1], save_resolution_scale*resolution[2]));
