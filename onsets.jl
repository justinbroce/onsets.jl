using  DSP
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

