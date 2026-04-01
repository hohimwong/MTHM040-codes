include("SJW_v5_LF.jl")
# script to run the SJW_v5_LF code

param = [-1000,1,1,1,1] # first element is dynamo number
bdy = [-1, 1]
δq = (bdy[2]-bdy[1])/ 300 
δt = 0.001
n = round(Int, ((bdy[2]-bdy[1])/δq) + 1 ) # number of grid points
m = 5000 # time steps
bdy_condi = [0,0,0,0]
#A0 = [(cos(i)) for i in LinRange(0,pi,n)] 
#B0 = [(cos(i)) for i in LinRange(0,pi,n)]
#A0 = [-i for i in LinRange(bdy[1],bdy[2],n)]
#B0 = [-i for i in LinRange(bdy[1],bdy[2],n)]
#A0 = [i^2 for i in LinRange(bdy[1],bdy[2],n)]
#B0 = [i^2 for i in LinRange(bdy[1],bdy[2],n)]
#A0 = [-i for i in LinRange(bdy[1],bdy[2],n)]
#B0 = [i^2 for i in LinRange(bdy[1],bdy[2],n)] 
A0 = [i^2 - i for i in LinRange(bdy[1],bdy[2],n)]
B0 = [i^2 - i for i in LinRange(bdy[1],bdy[2],n)]


results = SJW_v5_LF(param, δq, δt, bdy, m, A0, B0, bdy_condi, a1=0.1,a2=0.1,fq=1,EqA_sat_coeff=0,EqB_sat_coeff=1) # can adjust a1 a2 fq here for time dep.
A_results = copy(results[1])
B_results = copy(results[2])
transient_steps = 0 
# putting solution into matrices for contour plotting
B_results_matrix = zeros(Float64, m+1 - transient_steps, n)
for i in 1: (m+1 - transient_steps)
    B_results_matrix[i,:] = B_results[transient_steps + i]
end
A_results_matrix = zeros(Float64,m+1 - transient_steps,n)
for i in 1: (m+1-transient_steps)
    A_results_matrix[i,:] = A_results[transient_steps + i]
end

# plotting
fig1 = Figure()
fig2 = Figure()
ax1 = Axis(fig1[1,1],xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ],
ytickformat=yvalues ->["$(round(bdy[1]+(value-1)*δq, sigdigits=3))" for value in yvalues ],
xlabel = "t (δt = $(δt), $(m) steps)", ylabel = "q (q = cos(θ), δq = $(round(δq,sigdigits=4)))",
title = "Contour of B  (t ∈ [$(transient_steps*δt), $(m*δt)], q ∈ [$(bdy[1]), $(bdy[2])]), D = $(param[1]))")
contour!(ax1, B_results_matrix, levels = 55 , colormap = :hsv) 
ax2 = Axis(fig2[1,1],xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=3))" for value in values ], 
xlabel = "t (δt = $(δt), $(m) steps)", ylabel = "q (q = cos(θ), δq = $(round(δq, sigdigits=4)))",
ytickformat=yvalues ->["$(round(bdy[1]+(value-1)*δq, sigdigits=3))" for value in yvalues ], 
title = "Contour of A  (t ∈ [$(transient_steps*δt), $(m*δt)], q ∈ [$(bdy[1]), $(bdy[2])], D = $(param[1]))")
contour!(ax2, A_results_matrix, levels = 30, colormap=:hsv)
Colorbar(fig1[1,2], colormap =:hsv, limits = (findmin(B_results_matrix)[1], findmax(B_results_matrix)[1]), label = "value of B")
Colorbar(fig2[1,2], colormap=:hsv, limits = (findmin(A_results_matrix)[1], findmax(A_results_matrix)[1]), label = "value of A" )

fig_em = Figure()
ax_em = Axis(fig_em[1,1], title="Plot of magnetic energy over time", xlabel="time", ylabel="E_m",
xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ])
E_m = zeros(m+1-transient_steps)

for l in 1:length(E_m)[1] # loop to compute E_m at each time
    E_m[l] = sum((B_results[l+transient_steps].^2 + A_results[l+transient_steps].^2)*δq)*0.5
end

lines!(ax_em, E_m)

display(fig_em)
fig1

