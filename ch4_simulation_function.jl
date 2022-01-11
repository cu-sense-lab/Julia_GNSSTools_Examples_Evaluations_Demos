using Base.Threads

function simulate_signal(sigtype, f_s, t_length, prns, cn0s, f_ds, fd_rates, 
                         phis, code_start_idxs, h_parms; 
                         nADC=8, f_if=0, Tsys=535)
    @assert (length(prns) == length(cn0s) == length(f_ds) == length(fd_rates) ==
             length(phis) == length(code_start_idxs))
    if nADC == 8
        dtype = Complex{Int8}
    elseif nADC == 16
        dtype = Complex{Int16}
    elseif nADC == 32
        dtype = Complex{Int32}
    elseif nADC == 64
        dtype = Complex{Int64}
    else
        error("Invalid nADC value. Supported values are 8, 16, 32, and 64.")
    end
    sig_freq = sigtype.sig_freq
    sample_num = floor(Int, f_s*t_length)
    iq_values = zeros(Complex{Float64}, sample_num)
    signal = definesignal(sigtype, f_s, t_length; Tsys=Tsys, receiver_h_parms=h_parms)
    M = length(prns)
    # Add signal carriers
    for i in 1:M
        prn = prns[i]
        cn0 = cn0s[i]
        f_d = f_ds[i]
        fd_rate = fd_rates[i]
        phi = phis[i]
        code_start_idx = code_start_idxs[i]
        definesignal!(signal; prn=prn, f_if=f_if, f_d=f_d, fd_rate=fd_rate,
                      Tsys=Tsys, CN0=cn0, phi=phi, nADC=nADC, include_carrier=true,
                      include_thermal_noise=false, include_phase_noise=false,
                      include_adc=false, code_start_idx=code_start_idx)
        generatesignal!(signal)
        @threads for j in 1:sample_num
            iq_values[j] = iq_values[j] + signal.data[j]
        end
    end 
    B = maximum(skipmissing([sigtype.B_I, sigtype.B_Q]))
    σₙ = sqrt(GNSSTools.k*B*Tsys)
    # Add noise sources
    @threads for i in 1:sample_num
        iq_values[i] = iq_values[i]*cis(real(signal.phase_noise[i])) + σₙ*signal.thermal_noise[i]
    end
    # Quantize signal
    adc_scale = 2^(nADC-1)
    sigmax = sqrt(maximum(abs2, iq_values))
    @threads for i in 1:sample_num
        iq_values[i] = round(iq_values[i]*adc_scale/sigmax)
        if real(iq_values[i]) == adc_scale
            iq_values[i] = iq_values[i] - 1
        end
        if imag(iq_values[i]) == adc_scale
            iq_values[i] = iq_values[i] - 1im
        end
    end
    # Export GNSSTools.GNSSData struct
    data_type = string("sc", "$nADC")
    start_data_idx = 0
    data_start_time = missing
    site_loc_lla = missing
    total_data_length = t_length
    start_t = 0
    data = GNSSTools.GNSSData("", f_s, f_if, t_length, start_data_idx, 
                              dtype.(iq_values), data_type, data_start_time, 
                              site_loc_lla, sample_num, total_data_length, nADC, 
                              start_t)
    return data
end


function get_signal_truth_parms(t_length, f_s, code_start_idx, f_d, fd_rate, phi0, cn0,
                                Tsys, sigtype, channel, T; dt=0.001,
                                B=maximum(skipmissing([sigtype.B_I, sigtype.B_Q])))
    t = calctvector(floor(Int, t_length/dt), 1/dt)
    sig_freq = sigtype.sig_freq
    if channel == "Q"
        codes = sigtype.Q_codes
    elseif channel == "I"
        codes = sigtype.I_codes
    else
        error("Invalid channel value. Only 'I' and 'Q' allowed.")
    end
    code_length = codes.code_lengths[1]
    chipping_rate = codes.chipping_rates[1]
    phi = phi0 .+ 2π*(f_d .* t .+ 0.5 .* fd_rate .* t.^2)
    A = sqrt(2*GNSSTools.k*Tsys)*10^(cn0/20)
    σ = sqrt(GNSSTools.k*B*Tsys)
    PS = floor(Int, f_s*T)*A^2
    PN = σ^2
    snr = calc_snr(PS,PN)  
    doppler = f_d .+ fd_rate.*t
    doppler_rate = fill(fd_rate, length(t))
    f_code_d, f_code_dd = GNSSTools.calc_doppler_code_rate(chipping_rate, 
                                                           sig_freq, f_d, fd_rate)
    n0 = calcinitcodephase(code_length, f_code_d, f_code_dd, f_s,
                           code_start_idx)
    code_phase = n0 .+ f_code_d .* t .+ 0.5 .* f_code_dd .* t.^2
    return (t, doppler, doppler_rate, phi, code_phase, snr)
end