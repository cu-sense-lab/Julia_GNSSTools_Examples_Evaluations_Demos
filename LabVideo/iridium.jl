using GNSSTools
using SatelliteToolbox
using ProgressMeter
using Makie, FileIO, Colors, GeometryBasics, GLMakie
using AbstractPlotting, AbstractPlotting.MakieLayout


"""
Constellation info from 
[EOPortal](https://earth.esa.int/web/eoportal/satellite-missions/i/iridium-next)
"""


GM = 3.986004418e14  # m³s⁻²
eop = get_eop();
scale_a = Rₑ  + 1200e3  # meters


a = Rₑ+780*1000;  # meters
plane_num = 6;
sat_per_plane = 11;
incl = 86.4;
ΔΩ = 360/plane_num/2;
Δf_per_plane = 360/sat_per_plane/2;  # degrees
orbital_period = 2π*sqrt(a^3 / GM)  # seconds


Δt = 1 # seconds
t_range = Array(0:Δt:orbital_period);


constellation = define_constellation(a, plane_num, 
                                     sat_per_plane, incl, 
                                     t_range; obs_lla=missing,
                                     print_steps=true,
                                     eop=eop,
                                     ΔΩ=ΔΩ,
                                     show_plot=false,
                                     Δf_per_plane=Δf_per_plane);

# Earth sphere
u = Array(range(0, 2π, length=200)) .+ (180*pi/180);
v = Array(range(0, 1π, length=200));
u, v = meshgrid(u, v);
E_x = Rₑ.*cos.(u).*sin.(v);
E_y = Rₑ.*sin.(u).*sin.(v);
E_z = Rₑ.*cos.(v);


# Earth imagery
# earth = load(download("https://svs.gsfc.nasa.gov/vis/a000000/a002900/a002915/bluemarble-2048.png"));
earth = load("bluemarble-2048.png");


# View angle
ϕ = deg2rad(150)
θ = deg2rad(70)
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
        r_ecef = constellation.satellites[j].r_ecef[i,:]
        r_eci = kepler_to_sv(constellation.satellites[j].orbit[i]).r
        sat_pos[j,:] = r_eci
        if (constellation.satellites[j].id%sat_per_plane) == 1
            r_ecis = zeros(size(constellation.satellites[j].r_ecef))
            for k in 1:size(r_ecis)[1]
                r_ecis[k,:] =  kepler_to_sv(constellation.satellites[j].orbit[k]).r
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
                show_axis=false, lightposition=lightpos, 
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
update_cam!(scene, camera_pos, 0)
scene.center=false; save("iridium.png", scene, resolution=(save_resolution_scale*resolution[1], save_resolution_scale*resolution[2]));
