# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

include("../utils/read_nmvr.jl")

data_1 = read_nmvr("result.nmvr")
data_2 = read_nmvr("result.nmvr")

pressure_1 = data_1["cpx_pres"]
pressure_2 = data_2["cpx_pres"]

freq_1 = data_1["freq_stp"]
freq_2 = data_2["freq_stp"]

spl_1 = 20*log10.(abs.(pressure_1)/sqrt(2)/20e-6);
spl_2 = 20*log10.(abs.(pressure_2)/sqrt(2)/20e-6);

using Plots

Plots.plot(
    xaxis = :log10,
    xlim = [20, 100],
    # ylim = [0, 3],
    xticks = (20: 10 : 100, string.(20: 10 : 100)),
    xlabel = "Frequency (Hz)",
    ylabel = "SPL (dB)",
    title = "Frequency response (dB)",
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
    vec(freq_1), vec(spl_1[370,:]),
    label = "1"
)
Plots.plot!(
    vec(freq_2), vec(spl_2[370,:]),
    label = "2",
    linestyle = :dash
)
