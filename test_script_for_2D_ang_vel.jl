include("2D_ang_vel.jl")

# simple code to plot the spatially dependent u
x_bdy = [0,pi]
y_bdy = [0,0.5]
δx = (x_bdy[2] - x_bdy[1])/100
δy = (y_bdy[2] - y_bdy[1])/100
m_pts = round(Int, ((y_bdy[2] - y_bdy[1])/δy) + 1 )
n_pts = round(Int, ((x_bdy[2] - x_bdy[1])/δx) + 1 )

omega = TwoD_diff_rot(m_pts,n_pts,x_bdy,y_bdy,1)
f = Figure()
ax = Axis(f[1,1],title="u_space",xlabel="x", ylabel="z",
xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=3))" for value in values ],
ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=3))" for value in yvalues ],
aspect = 2.5)
display(omega)
surface!(ax,omega',colormap=:hsv)
Colorbar(f[1,2],colormap=:hsv,limits=(findmin(omega)[1],findmax(omega)[1]))
f