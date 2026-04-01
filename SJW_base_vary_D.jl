include("SJW_v5_LF.jl")
# simple script to run the SJW_v5_LF code for differnt values of D, and also plot the results
bdy = [-1, 1]
δq = (bdy[2]-bdy[1])/ 300 
δt = 0.001
n = round(Int, ((bdy[2]-bdy[1])/δq) + 1 ) # grid points
m = 5000 # time steps
bdy_condi = [0,0,0,0]

A0 = [i^2 - i for i in LinRange(bdy[1],bdy[2],n)]
B0 = [i^2 - i for i in LinRange(bdy[1],bdy[2],n)]

D_range = LinRange(0,5000,51)

for j in 1:length(D_range)
    param = [-D_range[j],1,1,1,1]

    results = SJW_v5_LF(param, δq, δt, bdy, m, A0, B0, bdy_condi, a1=0.5,a2=0.5,fq=0.5,EqA_sat_coeff=0,EqB_sat_coeff=1) # starting from amp=0.7 shows tailpiece
    A_results = copy(results[1])
    B_results = copy(results[2])
    transient_steps = 0
    B_results_matrix = zeros(Float64, m+1 - transient_steps, n)
    for i in 1: (m+1 - transient_steps)
        B_results_matrix[i,:] = B_results[transient_steps + i]
    end
    # plotting
    fig1 = Figure()
    ax1 = Axis(fig1[1,1],xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ],
    ytickformat=yvalues ->["$(round(bdy[1]+(value-1)*δq, sigdigits=3))" for value in yvalues ],
    xlabel = "t (δt = $(δt), $(m) steps)", ylabel = "q (q = cos(θ), δq = $(round(δq,sigdigits=4)))",
    title = "Contour of B  (t ∈ [$(transient_steps*δt), $(m*δt)], q ∈ [$(bdy[1]), $(bdy[2])]), D = $(round(Int,param[1])))")
    contour!(ax1, B_results_matrix, levels = 55 , colormap = :hsv)
    Colorbar(fig1[1,2], colormap =:hsv, limits = (findmin(B_results_matrix)[1], findmax(B_results_matrix)[1]), label = "value of B")

    fig_em = Figure()
    ax_em = Axis(fig_em[1,1], title="Plot of magnetic energy over time", xlabel="time", ylabel="E_m",
    xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ])
    E_m = zeros(m+1-transient_steps)
    for l in 1:length(E_m)[1]
        E_m[l] = sum((B_results[l+transient_steps].^2 + A_results[l+transient_steps].^2)*δq)*0.5
    end
    lines!(ax_em, E_m)
    display(fig_em)
    display(fig1)



end
