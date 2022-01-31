using  DSP, LinearAlgebra, ImageFilterting
#convolution version of moving average filter
# x is an array


function moving_average(x::AbstractArray,M::Integer)
    
    M > 1 || error("window size must be larger that 1")
    conv(x,rect(M)/M)[M ÷ 2 + 1 : length(x) + M÷2]

end

#adapted from https://www.audiolabs-erlangen.de/resources/MIR/FMP/C6/C6S1_NoveltyEnergy.html
function compute_novelty_energy(x, Fs=1, N=2048, H=128, gamma=10.0, norm=true)
    w = hanning(N)
    Fs_feature = Fs / H
    energy_local = conv(x.^2, w.^2,)
    energy_local = energy_local[1:H:length(energy_local)]
    if gamma != nothing
        energy_local = log.(1 .+ gamma .* energy_local)
    end
    energy_local_diff = diff(energy_local)
    push!(energy_local_diff,0)
    novelty_energy = copy(energy_local_diff)
    novelty_energy[energy_local_diff .< 0] .= 0
    if norm
        max_value = maximum(novelty_energy)
        if max_value > 0
            novelty_energy = novelty_energy / max_value
        end
    end
    return novelty_energy, Fs_feature
end
function peak_picking_simple(x, threshold=nothing)
 
    peaks = []
    if threshold == nothing
        threshold = minimum(x) - 1
    end
    for i in 2:length(x)
        if x[i - 1] < x[i] && x[i] > x[i + 1]
            if x[i] >= threshold
                push!(peaks,i)
            end
        end
    end
    return peaks
end

function median_filter(x,median_len)
    isodd(median_len) ? mapwindow(median, x, median_len)  : mapwindow(median, x, median_len+1) 
    
end
function peak_picking_MSAF(x, median_len=16, offset_rel=0.05, sigma=4.0)
    offset = mean(x) * offset_rel
    x = gaussian_filter(x, sigma)
    threshold_local = median_filter(x, median_len) .+ offset
    peaks = []
    for i in 2:length(x)-1
        if x[i - 1] < x[i] && x[i] > x[i + 1]
            if x[i] > threshold_local[i]
                push!(peaks,i)
            end
        end
    end
    
    return peaks, x, threshold_local
end
       
function gauss(σ::Real, l::Int = 4*ceil(Int,σ)+1)
    #isodd(l) || throw(ArgumentError("length must be odd"))
    w = l>>1
    g = σ == 0 ? [exp(0/(2*oftype(σ, 1)^2))] : [exp(-x^2/(2*σ^2)) for x=-w:w]
    centered(g/sum(g))
end
function gaussian_filter(x,σ)
    return (conv(gauss(σ),x)[1:length(x)])
    #==conv is slow
    and i was unable to convolve with filters
    if speed is an issue,
    replace code with imfilter from image filtering
    ==#
end

=##
        
