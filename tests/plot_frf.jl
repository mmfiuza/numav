
function read_matrix_binary(filename)
    open(filename, "r") do io
        rows = read(io, UInt64)
        cols = read(io, UInt64)
        n = rows * cols
        data = read!(io, Vector{ComplexF64}(undef, n))
        return reshape(data, (rows, cols))
    end
end

pressure_1 = read_matrix_binary("pressure_detj_fix_o2.bin")
pressure_2 = read_matrix_binary("pressure_area2_o2.bin")

spl_1 = 20*log10.(abs.(pressure_1)/sqrt(2)/20e-6);
spl_2 = 20*log10.(abs.(pressure_2)/sqrt(2)/20e-6);
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
    vec(freq), vec(spl_1[2,:]),
    label = "1"
)
Plots.plot!(
    vec(freq), vec(spl_2[2,:]),
    label = "2",
    linestyle = :dash
)
