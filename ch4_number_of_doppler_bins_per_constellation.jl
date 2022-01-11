using GNSSTools
using JLD

directory = "/home/sjbilardi/Dropbox/Apps/Overleaf/MastersThesis/"
doppler_file_name = string(directory, "ch1_doppler_curves.jld")
doppler_raw_gps, doppler_rate_raw_gps, sig_freq = load(doppler_file_name, 
                                                       "doppler_raw_gps",
                                                       "doppler_rate_raw_gps",
                                                       "freq")
doppler_raw_iridium, doppler_rate_raw_iridium = load(doppler_file_name, 
                                                     "doppler_raw_iridium",
                                                     "doppler_rate_raw_iridium")
doppler_raw_starlink, doppler_rate_raw_starlink = load(doppler_file_name, 
                                                       "doppler_raw_starlink",
                                                       "doppler_rate_raw_starlink")
doppler_raw_oneweb, doppler_rate_raw_oneweb = load(doppler_file_name, 
                                                   "doppler_raw_oneweb",
                                                   "doppler_rate_raw_oneweb")


Ts = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3]
constellation_dopplers = [doppler_raw_gps, 
                          doppler_raw_iridium,
                          doppler_raw_starlink,
                          doppler_raw_oneweb]
constellation_names = ["GPS            ", 
                       "Iridium        ", 
                       "Starlink       ", 
                       "OneWeb         "]
header = "Constellation"
for i in 1:length(Ts)
    global header = string(header, "\tT=$(floor(Int, Ts[i]*1000))")
end
println(header)
println("----------------------------------------------------------------------")
for i in 1:length(constellation_names)
    row = constellation_names[i]
    dopplers = constellation_dopplers[i]
    fd_range = max(maximum(dopplers), abs(minimum(dopplers)))
    for j in 1:length(Ts)
        Δf = 1/Ts[j]
        M = 2*ceil(Int, fd_range/Δf) + 1
        row = string(row, "\t$(M)")
    end
    println(row)
end
println("----------------------------------------------------------------------")

println("\n\n\nConstellation\tMax f_d (kHz)")
println("---------------------------------")
for i in 1:length(constellation_names)
    row = constellation_names[i]
    dopplers = constellation_dopplers[i]
    fd_range = max(maximum(dopplers), abs(minimum(dopplers)))
    fd_range = round(fd_range/1e3, digits=1)
    row = string(row, "\t$(fd_range)kHz")
    println(row)
end
println("------------------------------")
