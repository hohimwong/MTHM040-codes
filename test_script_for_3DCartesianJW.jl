include("3DCartesianJW.jl")

param = [-1,1,0.02]
δx = pi/100
δy = pi/100
δz = 0.01/2
δt = 0.05
x_bdy = [0,pi]
y_bdy = [0,pi]
z_bdy = [0.0, 0.5]
l_pts = round(Int, ((x_bdy[2] - x_bdy[1])/δx) + 1 )
m_pts = round(Int, ((y_bdy[2] - y_bdy[1])/δy) + 1 )
n_pts = round(Int, ((z_bdy[2] - z_bdy[1])/δz) + 1 )
i = 60
x_bdy_condi = zeros(6)
y_bdy_condi = zeros(6)
z_bdy_condi = zeros(6)

A0 = zeros(l_pts,m_pts,n_pts)
B0 = zeros(l_pts,m_pts,n_pts)
C0 = zeros(l_pts,m_pts,n_pts)

# the general form of the initial conditions below uses the IC from Kusano et al (2012), which is linear and force free (see below for the source)
# note this may not be that suitable of an IC for our problem; no results using this IC was obtained/discussed in the report 

for l in 1:l_pts
    for m in 1:m_pts
        for n in 1:n_pts

            A0[l,m,n] = (2*pi)^-1 * cos(2*pi*(m-1)*δy)*exp(-(4*pi^2 - 1)^0.5 * (n-1)*δz)
            B0[l,m,n] = -((4*pi^2 - 1)^0.5)*(1/(2*pi))*cos(2*pi*(m-1)*δy)*exp(-(4*pi^2 - 1)^0.5 * (n-1)*δz)
            C0[l,m,n] = sin(2*pi*(m-1)*δy)*exp(-(4*pi^2 - 1)^0.5 * (n-1)*δz)

        end
    end
end


ThreeDimCartesianJW(param, δx, δy, δz, δt, x_bdy, y_bdy, z_bdy, i, A0, B0, C0, x_bdy_condi, y_bdy_condi, z_bdy_condi, iter_steps=60)


# Kusano, K., Bamba, Y., Yamamoto, T. T., Iida, Y., Toriumi, S., & Asai, A. (2012). Magnetic field structures triggering solar flares and coronal mass ejections. The Astrophysical Journal, 760(1), 31.