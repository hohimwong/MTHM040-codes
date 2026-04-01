include("2D_JW_alpha2omega_v3.jl")

# script to run the 2D_JW_alpha2omega_v3 code, this code corresponds to choosing α and u that match JW91 (ϵ = 0)
# do note that "y" here actually corresponds to "z" in the actual model

param = [-20000,1,1,1]      # first element corresponds to D
x_bdy = [0,pi]
y_bdy = [0,0.5]
δt = 0.001
δx = (x_bdy[2] - x_bdy[1])/100 
δy = (y_bdy[2] - y_bdy[1])/100
m_pts = round(Int, ((y_bdy[2] - y_bdy[1])/δy) + 1 )
n_pts = round(Int, ((x_bdy[2] - x_bdy[1])/δx) + 1 ) # grid points in both spatial directions
i = 2000 # time steps
x_bdy_condi = [0,0,0,0]
y_bdy_condi = [0,0,0,0]
y_range = LinRange(y_bdy[1], y_bdy[2],m_pts)
A0 = zeros(m_pts,n_pts)
B0 = zeros(m_pts,n_pts)

for l in 1:m_pts
    A0[l,:] = [(sin(k)+cos(k)) for k in LinRange(x_bdy[1],x_bdy[2],n_pts)]
    B0[l,:] = [y_range[l]*(sin(k)+cos(k)) for k in LinRange(x_bdy[1],x_bdy[2],n_pts)]
end
#for l in 1:m_pts
#    A0[l,:] = [(sin(k)) for k in LinRange(x_bdy[1],x_bdy[2],n_pts)]
#    B0[l,:] = [(cos(k)) for k in LinRange(x_bdy[1],x_bdy[2],n_pts)]
#end 
#for l in 1:m_pts
#    A0[l,:] = [(cos(k)) for k in LinRange(x_bdy[1],x_bdy[2],n_pts)]
#    B0[l,:] = [(sin(k)) for k in LinRange(x_bdy[1],x_bdy[2],n_pts)]
#end 


results = TwoDJW_alp2ome_v3(param, δx, δy, δt, x_bdy, y_bdy, i, A0, B0, x_bdy_condi, y_bdy_condi, κ=1, τ=0, λ=1, ux=1, update=0,ϵ=0.0,tol=10e-3)
A_results = results[1]
B_results = results[2]
transient_steps = 0
B_results_matrix_y = zeros(i+1-transient_steps, n_pts)
B_results_matrix_x = zeros(i+1-transient_steps, n_pts)
x_cross_sec = 40
y_cross_sec = 50 # which cross sections to show
for l in 1:i+1-transient_steps
    B_results_matrix_y[l,:] = B_results[l+transient_steps][y_cross_sec,:]
end

for l in 1:i+1-transient_steps
    B_results_matrix_x[l,:] = B_results[l+transient_steps][:,x_cross_sec]
end

A_results_matrix = zeros(i+1-transient_steps, n_pts)
for l in 1:i+1-transient_steps
    A_results_matrix[l,:] = A_results[l+transient_steps][y_cross_sec,:]
end

# plotting
f = Figure()
ax = Axis(f[1,1],xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ],
ytickformat=yvalues ->["$(round(x_bdy[1]+(value-1)*δx, sigdigits=3))" for value in yvalues ],
xlabel = "t", ylabel="x, z = $((y_cross_sec-1)*δy)",
title = "Contour of B  (t ∈ [$(transient_steps*δt), $(i*δt)], x ∈ [$(x_bdy[1]), $(round(x_bdy[2], sigdigits=4))], z ∈ [$(y_bdy[1]), $(y_bdy[2])]), D = $(param[1]))")
contour!(ax, B_results_matrix_y, colormap =:hsv, levels=30)
Colorbar(f[1,2], colormap =:hsv, limits = (findmin(B_results_matrix_y)[1], findmax(B_results_matrix_y)[1]), label = "value of B")

f2 = Figure()
ax2 = Axis(f2[1,1],xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ],
ytickformat=yvalues ->["$(round(x_bdy[1]+(value-1)*δx, sigdigits=3))" for value in yvalues ],
xlabel = "t", ylabel="x (constant z)",
title = "Contour of A  (t ∈ [$(transient_steps*δt), $(i*δt)], x ∈ [$(x_bdy[1]), $(round(x_bdy[2], sigdigits=4))], z ∈ [$(y_bdy[1]), $(y_bdy[2])]), param = [$(param[1]), $(param[2]), $(param[3])])")
contour!(ax2, A_results_matrix, colormap =:viridis, levels=30)
Colorbar(f2[1,2], colormap=:viridis, limits=(findmin(A_results_matrix)[1], findmax(A_results_matrix)[1]), label="value of A")

E_m = zeros(i+1-transient_steps)
for l in 1:length(E_m)[1]
E_m[l] = sum((B_results[l+transient_steps].^2 + A_results[l+transient_steps].^2)*(δx*δy))*0.5
end

f3 = Figure()
ax3 = Axis(f3[1,1], ylabel="E_m", xlabel="t",
xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ], title="plot of magnetic energy")
lines!(ax3, E_m)

f4 = Figure()
ax4 = Axis(f4[1,1],xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ],
ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=3))" for value in yvalues ],
xlabel = "t", ylabel="y (we hold x constant)",
title = "Contour of B  (t ∈ [$(transient_steps*δt), $(i*δt)], x ∈ [$(x_bdy[1]), $(round(x_bdy[2], sigdigits=4))], y ∈ [$(y_bdy[1]), $(y_bdy[2])]), param = [$(param[1]), $(param[2]), $(param[3])])")
contour!(ax4, B_results_matrix_x, colormap=:hsv, levels = 30)
Colorbar(f4[1,2], colormap=:hsv, limits=(findmin(B_results_matrix_x)[1], findmax(B_results_matrix_x)[1]))

display(f3)
display(f2)
display(f)
