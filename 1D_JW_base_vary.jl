include("1D_JW_v4_LF.jl")
# simple script to run the 1D_JW_v4_LF code for differnt values of D, and also plot the results
bdy = [0,pi]
δx = (bdy[2]-bdy[1])/ 300 
δt = 0.001
n = round(Int, ((bdy[2]-bdy[1])/δx) + 1 ) # grid points
m = 5000 # time steps
bdy_condi = [0,0,0,0]

A0 = [(sin(i) + cos(i)) for i in LinRange(0,pi,n)]
B0 = [(sin(i) + cos(i)) for i in LinRange(0,pi,n)]

D_range = LinRange(0,5000,51)

for j in 1:length(D_range)
    param = [-D_range[j], 1, 1, 1, 1, 0, 1] 

    results = OneD_JW_v4_LF(param, δx, δt, bdy, m, A0, B0)
    A_results = copy(results[1])
    B_results = copy(results[2])
    transient_steps = 0
    B_results_matrix = zeros(Float64, m+1 - transient_steps, n)
    for i in 1: (m+1 - transient_steps)
        B_results_matrix[i,:] = B_results[transient_steps + i]
    end
    # plotting
    fig1 = Figure()
    ax1 = Axis(fig1[1,1], xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ],
    ytickformat=yvalues ->["$(round(bdy[1]+(value-1)*δx, sigdigits=3))" for value in yvalues ],
    xlabel = "t", ylabel="x",
    title = "Contour of B  (t ∈ [$(transient_steps*δt), $(m*δt)], x ∈ [$(bdy[1]), $(round(bdy[2], sigdigits=4))]), D = $(param[1]))")
    contour!(ax1, B_results_matrix, levels = 55 , colormap = :hsv)
    Colorbar(fig1[1,2], colormap =:hsv, limits = (findmin(B_results_matrix)[1], findmax(B_results_matrix)[1]), label = "value of B")

    fig_em = Figure()
    ax_em = Axis(fig_em[1,1], title="Plot of magnetic energy over time", xlabel="t", ylabel="E_m",
    xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ])
    E_m = zeros(m+1-transient_steps)
    for l in 1:length(E_m)[1]
        E_m[l] = sum((B_results[l+transient_steps].^2 + A_results[l+transient_steps].^2)*δx)*0.5
    end
    lines!(ax_em, E_m)
    display(fig_em)
    display(fig1)



end
