include("3DJW_TrueSoln_compare.jl")
# script to run the 3DJW_TrueSoln_compare code

param = [-1,1,1]
δx = pi/100
δy = pi/100
δz = 0.01/2
δt = 0.025
x_bdy = [0,pi]
y_bdy = [0,pi]
z_bdy = [0.0, 0.5]
l_pts = round(Int, ((x_bdy[2] - x_bdy[1])/δx) + 1 )
m_pts = round(Int, ((y_bdy[2] - y_bdy[1])/δy) + 1 )
n_pts = round(Int, ((z_bdy[2] - z_bdy[1])/δz) + 1 ) # grid points
i = 10 # time steps
x_bdy_condi = zeros(6)
y_bdy_condi = zeros(6)
z_bdy_condi = zeros(6)

A0 = zeros(l_pts,m_pts,n_pts)
B0 = zeros(l_pts,m_pts,n_pts)
C0 = zeros(l_pts,m_pts,n_pts)

for l in 1:l_pts
    for m in 1:m_pts
        for n in 1:n_pts
            A0[l,m,n] = sin((l-1)*δx)*sin((m-1)*δy)*sin(2*pi*(n-1)*δz)
            B0[l,m,n] = sin((l-1)*δx)*sin((m-1)*δy)*sin(2*pi*(n-1)*δz)
            C0[l,m,n] = sin((l-1)*δx)*sin((m-1)*δy)*sin(2*pi*(n-1)*δz)
        end
    end
end

ThreeDimJWTrueSoln(param, δx, δy, δz, δt, x_bdy, y_bdy, z_bdy, i, A0, B0, C0, x_bdy_condi, y_bdy_condi, z_bdy_condi, iter_steps=50, line_step=150, vec_field="false")