
function read_matrix_binary(filename)
    open(filename, "r") do io
        rows = read(io, UInt64)
        cols = read(io, UInt64)
        n = rows * cols
        data = read!(io, Vector{ComplexF64}(undef, n))
        return reshape(data, (rows, cols))
    end
end

pressure = read_matrix_binary("pressure.bin")

spl = 20*log10.(abs.(pressure)/sqrt(2)/20e-6);
freq = LinRange(20, 100, 1000)

using Plots

Plots.plot(
    xaxis = :log10,
    xlim = [20, 100],
    ylim = [0, 200],
    xticks = (20: 10 : 100, string.(20: 10 : 100)),
    xlabel = "Frequency (Hz)",
    ylabel = "SPL (dB)",
    title = "Frequency response",
    grid = true,
    gridalpha = 0.25,
    legendposition = (0.11, 0.9),
    framestyle = :box,
    linestyle = :solid,
    titlefontsize = 14,
    legendfontsize = 10,
    xtickfontsize = 10,
    ytickfontsize = 10,
    xguidefontsize = 12,
    yguidefontsize = 12,
    size = (600, 400),
)
Plots.plot!(
    vec(freq), vec(spl[100,:]),
    label = "FRF"
)
