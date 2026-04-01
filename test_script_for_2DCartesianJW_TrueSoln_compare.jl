include("2DCartesianJW_TrueSol_compare.jl")
# script to run the code 2DCartesianJW_TrueSol_compare
# note the y's correspond to z
param = [-1,1,1] # dynamo number does not matter here since we are only dealing with the diffusion component
x_bdy = [0,pi]
y_bdy = [0,0.5]
δt = 0.001
δx = (x_bdy[2] - x_bdy[1])/50
δy = (y_bdy[2] - y_bdy[1])/50
m_pts = round(Int, ((y_bdy[2] - y_bdy[1])/δy) + 1 )
n_pts = round(Int, ((x_bdy[2] - x_bdy[1])/δx) + 1 ) # grid points in both directions
i = 300 # time steps
x_bdy_condi = [0,0,0,0]
y_bdy_condi = [0,0,0,0]
A0 = zeros(m_pts,n_pts)
B0 = zeros(m_pts,n_pts)
x_range = LinRange(x_bdy[1],x_bdy[2],n_pts)
y_range = LinRange(y_bdy[1],y_bdy[2],m_pts)

for l in 1:m_pts
    for s in 1:n_pts
        A0[l,s] = sin(x_range[s])*sin(2*pi*y_range[l]) # initial condition corresponding to the true solution in the code
    end
end
B0 = copy(A0)

results = TwoDimCartesianJW_Truesol_comp(param, δx, δy, δt, x_bdy, y_bdy, i, A0, B0, x_bdy_condi, y_bdy_condi, κ=1, τ=0, λ=1, iter_steps=50, source_term_switch=0, update=1)
A_results = results[1]
B_results = results[2]
error_list = results[3]
max_numer_soln_list = results[4]
max_true_soln_list = results[5]
transient_steps = 0
B_results_matrix_y = zeros(i+1-transient_steps, n_pts)
B_results_matrix_x = zeros(i+1-transient_steps, n_pts)
x_cross_sec = 40
y_cross_sec = 40 # can be changes to the desired index to show a differnt cross section
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

f5 = Figure()
ax5 = Axis(f5[1,1], xlabel="t", ylabel="log of Maximum error |Aₜᵣᵤₑ-A|",xtickformat=values -> ["$(round(δt*(value+(transient_steps-1)), sigdigits=4))" for value in values ])
lines!(ax5, log.(error_list), color=:"blue", label="Maximum error")


f5