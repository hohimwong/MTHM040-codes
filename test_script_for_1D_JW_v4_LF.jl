include("1D_JW_v4_LF.jl")

# script to run the 1D_JW_v4_LF code
param = [-1000, 1, 1, 1, 1, 0, 1]   # the first element is the dynamo number
bdy = [0,pi]
m = 50000 # number of timesteps
δx = pi/300
δt = 0.001
n = round(Int, ((bdy[2]-bdy[1])/δx) +1) # number of grid points
#A0 = [(sin(i)) for i in LinRange(0,pi,n)]
#B0 = [(sin(i)) for i in LinRange(0,pi,n)]
A0 = [(cos(i)) for i in LinRange(0,pi,n)]
B0 = [(sin(i)) for i in LinRange(0,pi,n)]
#A0 = [(sin(i) + cos(i)) for i in LinRange(0,pi,n)]
#B0 = [(sin(i) + cos(i)) for i in LinRange(0,pi,n)] 
results = OneD_JW_v4_LF(param, δx, δt, bdy, m, A0, B0)
A_results = copy(results[1])
B_results = copy(results[2])

transient_steps = 0
B_results_matrix = zeros(Float64, m+1 - transient_steps, n) # putting the results into matrices to plot the contours
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
ax1 = Axis(fig1[1,1], xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ],
ytickformat=yvalues ->["$(round(bdy[1]+(value-1)*δx, sigdigits=3))" for value in yvalues ],
xlabel = "t", ylabel="x",
title = "Contour of B  (t ∈ [$(transient_steps*δt), $(m*δt)], x ∈ [$(bdy[1]), $(round(bdy[2], sigdigits=4))]), D = $(param[1]))")

contour!(ax1, B_results_matrix, levels = 55 , colormap = :hsv)
ax2 = Axis(fig2[1,1], xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ],
ytickformat=yvalues ->["$(round(bdy[1]+(value-1)*δx, sigdigits=3))" for value in yvalues ],
xlabel = "t", ylabel="x",
title = "Contour of A  (t ∈ [$(transient_steps*δt), $(m*δt)], x ∈ [$(bdy[1]), $(round(bdy[2], sigdigits=4))]), D = $(param[1]))")

contour!(ax2, A_results_matrix, levels = 55, colormap=:hsv)
Colorbar(fig1[1,2], colormap =:hsv, limits = (findmin(B_results_matrix)[1], findmax(B_results_matrix)[1]), label = "value of B")
Colorbar(fig2[1,2], colormap=:hsv, limits = (findmin(A_results_matrix)[1], findmax(A_results_matrix)[1]), label = "value of A" )

fig_em = Figure()
ax_em = Axis(fig_em[1,1], title="Plot of magnetic energy over time", xlabel="t", ylabel="E_m",
xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ])
E_m = zeros(m+1-transient_steps)

for l in 1:length(E_m)[1] # loop to compute E_m at each time
    E_m[l] = sum((B_results[l+transient_steps].^2 + A_results[l+transient_steps].^2)*δx)*0.5
end

lines!(ax_em, E_m)
display(fig_em)
fig1
# can type the names of the figure to show it, e.g. typing "fig1" in the terminal shows the contours for B