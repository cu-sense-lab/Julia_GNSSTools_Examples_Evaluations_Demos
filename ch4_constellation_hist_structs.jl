using Statistics
using Random
using StatsBase
using LinearAlgebra


struct ConstellationHistogram{T1,T2,T3,T4,T5}
    nbins::Int
    dopplers::Array{Float64,1}
    doppler_rates::Array{Float64,1}
    histogram::T1
    weights::T2
    fd_rate::Vector{Float64}
    f_d::Vector{Float64}
    elevations::T3
    elev_hists::T4
    max_elevations::T5
end


function get_histogram(x; nbins)
    return normalize(fit(Histogram, x, nbins=nbins), mode=:probability)
end


function initiate_hist(dopplers, doppler_rates; nbins=1000)
    histogram = normalize(fit(Histogram, (doppler_rates, dopplers), 
                              nbins=nbins), mode=:probability)
    weights = collect(Iterators.flatten(histogram.weights))
    weights = Weights(weights)
    fd_rate, f_d = meshgrid(Array(histogram.edges[1][1:end-1]), 
                            Array(histogram.edges[2][1:end-1]))
    fd_rate = collect(Iterators.flatten(fd_rate))
    f_d = collect(Iterators.flatten(f_d))
    return ConstellationHistogram(nbins, dopplers, doppler_rates, 
                                  histogram, weights, fd_rate, f_d, 
                                  missing, missing, missing)
end


function initiate_hist(dopplers, doppler_rates, elevations, max_elevations; 
                       nbins=1000, elev_bins=100)
    histogram = normalize(fit(Histogram, (doppler_rates, dopplers), 
                              nbins=nbins), mode=:probability)
    histogram3d = normalize(fit(Histogram, (elevations, doppler_rates, dopplers), 
                               nbins=(elev_bins, nbins, nbins)), mode=:probability)
    histogram_dopplers =  normalize(fit(Histogram, (elevations, dopplers), 
                                        nbins=(elev_bins, nbins)), mode=:probability)
    histogram_doppler_rates =  normalize(fit(Histogram, (elevations, doppler_rates), 
                                        nbins=(elev_bins, nbins)), mode=:probability)
    elev_hists = [histogram3d, histogram_dopplers, histogram_doppler_rates]                          
    weights = collect(Iterators.flatten(histogram.weights))
    weights = Weights(weights)
    fd_rate, f_d = meshgrid(Array(histogram.edges[1][1:end-1]), 
                            Array(histogram.edges[2][1:end-1]))
    fd_rate = collect(Iterators.flatten(fd_rate))
    f_d = collect(Iterators.flatten(f_d))
    return ConstellationHistogram(nbins, dopplers, doppler_rates, 
                                  histogram, weights, fd_rate, f_d, 
                                  elevations, elev_hists, max_elevations)
end


function sample_distribution(histogram::ConstellationHistogram) 
    f_d = sample(histogram.f_d, histogram.weights)
    fd_rate = sample(histogram.fd_rate, histogram.weights)
    return (f_d, fd_rate)
end


function max_doppler_curve_per_elev(histogram::ConstellationHistogram, elevs)
    elevations = histogram.elevations
    
    doppler_bins = histogram.histogram.edges[1]
    dopplers = histogram.dopplers
    doppler_rates = histogram.doppler_rates
    N = size(histogram.elev_hists[2].weights)[1]
    doppler_medians = zeros(N)
    doppler_rate_medians = zeros(N)
    for i in 1:N
        doppler_median_idx = argmax(histogram.elev_hists[2].weights[i,:])
        doppler_rate_median_idx = argmax(histogram.elev_hists[3].weights[i,:])
        doppler_medians[i] = histogram.elev_hists[2].edges[2][doppler_median_idx]
        doppler_rate_medians[i] = histogram.elev_hists[3].edges[2][doppler_rate_median_idx]
    end
    return (doppler_medians, doppler_rate_medians)
end