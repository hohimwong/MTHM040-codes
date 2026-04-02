using LinearAlgebra
using Statistics
using CairoMakie
# slightly altered version of the 3D code, where we test the diffusion component
# the results do not fully match the expectations; there are likely errors in the code

function ADI_eqn_A_x!(A_old, A_new, B_old, C_old, l_pts, m_pts, n_pts, x_bdy_condi, a1, a2, a3, a4, a5, a6, a7)
    # using lower case a b c to denote the coefficients
    
    r = zeros(l_pts)

    M_A_x = zeros(n_pts, l_pts)
    M_A_x[1,1] = 1
    M_A_x[end,end] = 1
    for j in 2:n_pts-1
        M_A_x[j,j-1] = -a1
        M_A_x[j,j] = a4 - a5[j,j] + 1
        M_A_x[j,j+1] = -a1
    end
    LHS_M_inv = inv(M_A_x)

    for m in 2:m_pts-1 # sweep in x direction
        for n in 2:n_pts-1
            r[1] = x_bdy_condi[1]
            for l in 2:l_pts-1 
                r[l] = (-a5[l,n]+a2)*A_new[l,m+1,n] + (a2)*A_new[l,m-1,n] + (a3)*A_new[l,m,n+1] + (a3)*A_new[l,m,n-1] + 
                    (a1)*A_old[l+1,m,n] + (a5[l,n]-a4+1)*A_old[l,m,n] + (a1)*A_old[l-1,m,n] + (-a5[l,n]+a2)*A_old[l,m+1,n] +
                        (a2)*A_old[l,m-1,n] + (a3)*A_old[l,m,n+1] + (a3)*A_old[l,m,n-1] + (a6[l])*C_old[l,m+1,n] -
                            (a6[l])*C_old[l,m,n] - (a7[l])*B_old[l,m,n+1] + (a7[l])*B_old[l,m,n]
            end
            r[end] = x_bdy_condi[2]
            A_new[:,m,n] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return A_new

end

function ADI_eqn_A_y!(A_old, A_new, B_old, C_old, l_pts, m_pts, n_pts, y_bdy_condi, a1, a2, a3, a4, a5, a6, a7)

    r = zeros(m_pts)
    M_A_y = zeros(n_pts, m_pts)

    for l in 2:l_pts-1 # y direction
        M_A_y[1,1] = 1
        M_A_y[end,end] = 1
        for j in 2:n_pts-1
            M_A_y[j,j-1] = -a2
            M_A_y[j,j] = a4 - a5[l,j] + 1
            M_A_y[j,j+1] = a5[l,j+1] - a2
        end
        LHS_M_inv = inv(M_A_y)

        for n in 2:n_pts-1
            r[1] = y_bdy_condi[1]
            for m in 2:m_pts-1 
                r[m] = (a1)*A_new[l+1,m,n] + (a1)*A_new[l-1,m,n] + (a3)*A_new[l,m,n+1] + (a3)*A_new[l,m,n-1] + 
                    (a1)*A_old[l+1,m,n] + (a5[l,n]-a4+1)*A_old[l,m,n] + (a1)*A_old[l-1,m,n] + (-a5[l,n]+a2)*A_old[l,m+1,n] +
                        (a2)*A_old[l,m-1,n] + (a3)*A_old[l,m,n+1] + (a3)*A_old[l,m,n-1] + (a6[l])*C_old[l,m+1,n] -
                            (a6[l])*C_old[l,m,n] - (a7[l])*B_old[l,m,n+1] + (a7[l])*B_old[l,m,n]
            end
            r[end] = y_bdy_condi[2]
            A_new[l,:,n] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return A_new

end

function ADI_eqn_A_z!(A_old, A_new, B_old, C_old, l_pts, m_pts, n_pts, z_bdy_condi, a1, a2, a3, a4, a5, a6, a7)
    
    r = zeros(n_pts)
    M_A_z = zeros(l_pts, n_pts)
    M_A_z[1,1] = 1
    M_A_z[end,end] = 1
    for j in 2:n_pts-1
        M_A_z[j,j-1] = -a3
        M_A_z[j,j] = a4 - a5[j,j] + 1
        M_A_z[j,j+1] = -a3
    end
    LHS_M_inv = inv(M_A_z)

    for m in 2:m_pts-1  # z direction
        for l in 2:l_pts-1
            r[1] = z_bdy_condi[1]
            for n in 2:n_pts-1 
                r[n] = (-a5[l,n]+a2)*A_new[l,m+1,n] + (a2)*A_new[l,m-1,n] + (a1)*A_new[l+1,m,n] + (a1)*A_new[l-1,m,n] + 
                    (a1)*A_old[l+1,m,n] + (a5[l,n]-a4+1)*A_old[l,m,n] + (a1)*A_old[l-1,m,n] + (-a5[l,n]+a2)*A_old[l,m+1,n] +
                        (a2)*A_old[l,m-1,n] + (a3)*A_old[l,m,n+1] + (a3)*A_old[l,m,n-1] + (a6[l])*C_old[l,m+1,n] -
                            (a6[l])*C_old[l,m,n] - (a7[l])*B_old[l,m,n+1] + (a7[l])*B_old[l,m,n]
            end
            r[end] = z_bdy_condi[2]
            A_new[l,m,:] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return A_new

end

function ADI_eqn_B_x!(A_old, A_new, B_old, B_new, C_old, l_pts, m_pts, n_pts, x_bdy_condi, b1, b2, b3, b4, b5, b6, b7, b8, b9)
    
    r = zeros(l_pts)
    M_B_x = zeros(n_pts, l_pts)
    M_B_x[1,1] = 1
    M_B_x[end,end] = 1
    for j in 2:n_pts-1
        M_B_x[j,j-1] = -b1
        M_B_x[j,j] = 1 + b4 - b5[j,j]
        M_B_x[j,j+1] = -b1
    end
    LHS_M_inv = inv(M_B_x)

    for m in 2:m_pts-1 # x direction
        for n in 2:n_pts-1
            r[1] = x_bdy_condi[3]
            for l in 2:l_pts-1 
                r[l] = (-b5[l,n]+b2)*B_new[l,m+1,n] + (b2)*B_new[l,m-1,n] + (b3)*B_new[l,m,n+1] + (b3)*B_new[l,m,n-1] + 
                    (b1)*B_old[l+1,m,n] + (b5[l,n]-b4+1)*B_old[l,m,n] + (b1)*B_old[l-1,m,n] + (-b5[l,n]+b2)*B_old[l,m+1,n] + 
                        (b2)*B_old[l,m-1,n] + (b3)*B_old[l,m,n+1] + (b3)*B_old[l,m,n-1] + (b6[l])*C_old[l,m,n] + 
                            (b7[l])*A_new[l,m,n+1] - (b7[l])*A_new[l,m,n] + (b7[l])*A_old[l,m,n+1] - (b7[l])*A_old[l,m,n] - 
                                (b8[l+1])*C_old[l+1,m,n] + (b8[l])*C_old[l,m,n] - (b9[l])*C_old[l,m,n]
            end
            r[end] = x_bdy_condi[4]
            B_new[:,m,n] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return B_new

end

function ADI_eqn_B_y!(A_old, A_new, B_old, B_new, C_old, l_pts, m_pts, n_pts, y_bdy_condi, b1, b2, b3, b4, b5, b6, b7, b8, b9)
    
    r = zeros(m_pts)
    M_B_y = zeros(n_pts, m_pts)

    for l in 2:l_pts-1 # y direction
        M_B_y[1,1] = 1
        M_B_y[end,end] = 1
        for j in 2:n_pts-1
            M_B_y[j,j-1] = -b2
            M_B_y[j,j] = b4 - b5[l,j] + 1
            M_B_y[j,j+1] = b5[l,j+1] - b2
        end
        LHS_M_inv = inv(M_B_y)

        for n in 2:n_pts-1
            r[1] = y_bdy_condi[3]
            for m in 2:m_pts-1 
                r[m] = (b1)*B_new[l+1,m,n] + (b1)*B_new[l-1,m,n] + (b3)*B_new[l,m,n+1] + (b3)*B_new[l,m,n-1] + 
                    (b1)*B_old[l+1,m,n] + (b5[l,n]-b4+1)*B_old[l,m,n] + (b1)*B_old[l-1,m,n] + (-b5[l,n]+b2)*B_old[l,m+1,n] + 
                        (b2)*B_old[l,m-1,n] + (b3)*B_old[l,m,n+1] + (b3)*B_old[l,m,n-1] + (b6[l])*C_old[l,m,n] + 
                            (b7[l])*A_new[l,m,n+1] - (b7[l])*A_new[l,m,n] + (b7[l])*A_old[l,m,n+1] - (b7[l])*A_old[l,m,n] - 
                                (b8[l+1])*C_old[l+1,m,n] + (b8[l])*C_old[l,m,n] - (b9[l])*C_old[l,m,n]
            end
            r[end] = y_bdy_condi[4]
            B_new[l,:,n] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return B_new

end

function ADI_eqn_B_z!(A_old, A_new, B_old, B_new, C_old, l_pts, m_pts, n_pts, z_bdy_condi, b1, b2, b3, b4, b5, b6, b7, b8, b9)
    
    r = zeros(n_pts)
    M_B_z = zeros(l_pts, n_pts)
    M_B_z[1,1] = 1
    M_B_z[end,end] = 1
    for j in 2:n_pts-1
        M_B_z[j,j-1] = -b3
        M_B_z[j,j] = b4 - b5[j,j] + 1
        M_B_z[j,j+1] = -b3
    end
    LHS_M_inv = inv(M_B_z)
    
    for m in 2:m_pts-1 # z direction
        for l in 2:l_pts-1
            r[1] = z_bdy_condi[3]
            for n in 2:n_pts-1 
                r[n] = (-b5[l,n]+b2)*B_new[l,m+1,n] + (b2)*B_new[l,m-1,n] + (b1)*B_new[l+1,m,n] + (b1)*B_new[l+1,m,n] + 
                    (b1)*B_old[l+1,m,n] + (b5[l,n]-b4+1)*B_old[l,m,n] + (b1)*B_old[l-1,m,n] + (-b5[l,n]+b2)*B_old[l,m+1,n] + 
                        (b2)*B_old[l,m-1,n] + (b3)*B_old[l,m,n+1] + (b3)*B_old[l,m,n-1] + (b6[l])*C_old[l,m,n] + 
                            (b7[l])*A_new[l,m,n+1] - (b7[l])*A_new[l,m,n] + (b7[l])*A_old[l,m,n+1] - (b7[l])*A_old[l,m,n] - 
                                (b8[l+1])*C_old[l+1,m,n] + (b8[l])*C_old[l,m,n] - (b9[l])*C_old[l,m,n]
            end
            r[end] = z_bdy_condi[4]
            B_new[l,m,:] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return B_new

end

function ADI_eqn_C_x!(A_old, A_new, B_old, B_new, C_old, C_new, l_pts, m_pts, n_pts, x_bdy_condi, c1, c2, c3, c4, c5, c6, c7, c8)
    
    r = zeros(l_pts)
    M_C_x = zeros(n_pts, l_pts)
    M_C_x[1,1] = 1
    M_C_x[end,end] = 1
    for j in 2:n_pts-1
        M_C_x[j,j-1] = -c1
        M_C_x[j,j] = c4 - c5[j,j] + 1
        M_C_x[j,j+1] = -c1
    end
    LHS_M_inv = inv(M_C_x)

    for m in 2:m_pts-1 # x direction
        for n in 2:n_pts-1
            r[1] = x_bdy_condi[5]
            for l in 2:l_pts-1 
                r[l] = (-c5[l,n]+c2)*C_new[l,m+1,n] + (c2)*C_new[l,m-1,n] + (c3)*C_new[l,m,n+1] + (c3)*C_new[l,m,n-1] + 
                    (c1)*C_old[l+1,m,n] + (c5[l,n]-c4+1)*C_old[l,m,n] + (c1)*C_old[l-1,m,n] + (-c5[l,n]+c2)*C_old[l,m+1,n] + (c2)*C_old[l,m-1,n] +
                        (c3)*C_old[l,m,n+1] + (c3)*C_old[l,m,n-1] + (c6[l])*B_new[l,m,n] + (c6[l])*B_old[l,m,n] + (c7[l+1])*B_new[l+1,m,n] - 
                            (c7[l])*B_new[l,m,n] + (c7[l+1])*B_old[l+1,m,n] - (c7[l])*B_old[l,m,n] - (c8[l])*A_new[l,m+1,n] + (c8[l])*A_new[l,m,n] -
                                (c8[l])*A_old[l,m+1,n] + (c8[l])*A_old[l,m,n]
            end
            r[end] = x_bdy_condi[6]
            C_new[:,m,n] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return C_new

end

function ADI_eqn_C_y!(A_old, A_new, B_old, B_new, C_old, C_new, l_pts, m_pts, n_pts, y_bdy_condi, c1, c2, c3, c4, c5, c6, c7, c8)
    
    r = zeros(m_pts)
    M_C_y = zeros(n_pts, m_pts)

    for l in 2:l_pts-1 # y direction
        M_C_y[1,1] = 1
        M_C_y[end,end] = 1
        for j in 2:n_pts-1
            M_C_y[j,j-1] = -c2
            M_C_y[j,j] = c4 - c5[l,j] + 1
            M_C_y[j,j+1] = c5[l,j+1] - c2
        end
        LHS_M_inv = inv(M_C_y)

        for n in 2:n_pts-1
            r[1] = y_bdy_condi[5]
            for m in 2:m_pts-1 
                r[m] = (c1)*C_new[l+1,m,n] + (c1)*C_new[l-1,m,n] + (c3)*C_new[l,m,n+1] + (c3)*C_new[l,m,n-1] + 
                    (c1)*C_old[l+1,m,n] + (c5[l,n]-c4+1)*C_old[l,m,n] + (c1)*C_old[l-1,m,n] + (-c5[l,n]+c2)*C_old[l,m+1,n] + (c2)*C_old[l,m-1,n] +
                        (c3)*C_old[l,m,n+1] + (c3)*C_old[l,m,n-1] + (c6[l])*B_new[l,m,n] + (c6[l])*B_old[l,m,n] + (c7[l+1])*B_new[l+1,m,n] - 
                            (c7[l])*B_new[l,m,n] + (c7[l+1])*B_old[l+1,m,n] - (c7[l])*B_old[l,m,n] - (c8[l])*A_new[l,m+1,n] + (c8[l])*A_new[l,m,n] -
                                (c8[l])*A_old[l,m+1,n] + (c8[l])*A_old[l,m,n]
            end
            r[end] = y_bdy_condi[6]
            C_new[l,:,n] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return C_new

end

function ADI_eqn_C_z!(A_old, A_new, B_old, B_new, C_old, C_new, l_pts, m_pts, n_pts, z_bdy_condi, c1, c2, c3, c4, c5, c6, c7, c8)
    
    r = zeros(n_pts)
    M_C_z = zeros(l_pts, n_pts)
    M_C_z[1,1] = 1
    M_C_z[end,end] = 1
    for j in 2:n_pts-1
        M_C_z[j,j-1] = -c3
        M_C_z[j,j] = c4 - c5[j,j] + 1
        M_C_z[j,j+1] = -c3
    end
    LHS_M_inv = inv(M_C_z)
    for m in 2:m_pts-1 # z direction
        for l in 2:l_pts-1
            r[1] = z_bdy_condi[5]
            for n in 2:n_pts-1 
                r[n] = (-c5[l,n]+c2)*C_new[l,m+1,n] + (c2)*C_new[l,m-1,n] + (c1)*C_new[l+1,m,n] + (c1)*C_new[l-1,m,n] + 
                    (c1)*C_old[l+1,m,n] + (c5[l,n]-c4+1)*C_old[l,m,n] + (c1)*C_old[l-1,m,n] + (-c5[l,n]+c2)*C_old[l,m+1,n] + (c2)*C_old[l,m-1,n] +
                        (c3)*C_old[l,m,n+1] + (c3)*C_old[l,m,n-1] + (c6[l])*B_new[l,m,n] + (c6[l])*B_old[l,m,n] + (c7[l+1])*B_new[l+1,m,n] - 
                            (c7[l])*B_new[l,m,n] + (c7[l+1])*B_old[l+1,m,n] - (c7[l])*B_old[l,m,n] - (c8[l])*A_new[l,m+1,n] + (c8[l])*A_new[l,m,n] -
                                (c8[l])*A_old[l,m+1,n] + (c8[l])*A_old[l,m,n]
            end
            r[end] = z_bdy_condi[6]
            C_new[l,m,:] = copy(LHS_M_inv * r) # matrix depends on the LHS coefficients
        end
    end

    return C_new

end

function ThreeDimJWTrueSoln(param, δx, δy, δz, δt, x_bdy, y_bdy, z_bdy, i, A0, B0, C0, x_bdy_condi, y_bdy_condi, z_bdy_condi; iter_steps=20, line_step=30, vec_field="true", contour="true",source="false")

    α0 = copy(param[1])
    ω0 = copy(param[2])
    η = copy(param[3])
    l_pts = round(Int, ((x_bdy[2] - x_bdy[1])/δx) + 1)
    display(l_pts)
    m_pts = round(Int, ((y_bdy[2] - y_bdy[1])/δy) + 1)
    display(m_pts)
    n_pts = round(Int, ((z_bdy[2] - z_bdy[1])/δz) + 1)
    display(n_pts)

    x_range = LinRange(x_bdy[1], x_bdy[2], l_pts)
    y_range = LinRange(y_bdy[1], y_bdy[2], m_pts)
    z_range = LinRange(z_bdy[1], z_bdy[2], n_pts)

    ψ = 0
    if source == "false"
        ψ = 0
    elseif source == "true"
        ψ = 1
    end

    # First we create the vectors and matrices for α ω and u 

    α = α0 * cos.(x_range)
    αₓ = α0 * sin.(x_range)
    u = zeros(l_pts,n_pts)
    for j in 1:l_pts
        for k in 1:n_pts
            u[j,k] = ω0 * z_range[k] #* sin(x_range[j])
        end
    end
    ω = ω0 *ones(l_pts)#* sin.(x_range)

    # we define the coefficients that will be used later in the ADI equations
    # The k coefficients are the same for all ADI equations, and the ones denoted by a b c corresponds to the A B C equations

    # k1 to k4 are diffusion terms, else are source terms, which will be turned off
    k1 = (η*δt)/(2*(δx^2))
    k2 = (η*δt)/(2*(δy^2))
    k3 = (η*δt)/(2*(δz^2))
    k4 = (η*δt)/(δx^2) + (η*δt)/(δy^2) + (η*δt)/(δz^2)
    k5 = ψ * u * (δt)/(2*δy)        # this is a matrix

    # The coefficients for A

    a6 = ψ * α * (δt)/(δy)          # this is a vector
    a7 = ψ * α * (δt)/(δz)          # this is a vector

    # The coefficients for B 

    b6 = ψ * δt * ω                 # this is a vector
    b7 = ψ * α * (δt)/(2*δz)        # this is a vector
    b8 = ψ * α * (δt)/(δx)          # this is a vector
    b9 = ψ * δt * αₓ                # this is a vector

    # The coefficients for C

    c6 = ψ * αₓ * (δt)/(2)
    c7 = ψ * α * (δt)/(2*δx)
    c8 = ψ * α * (δt)/(2*δy)

    # Then make sure the boundary conditions are applied for the initial state

    A0[1,:,:] .= x_bdy_condi[1]
    A0[end,:,:] .= x_bdy_condi[2]
    A0[:,1,:] .= y_bdy_condi[1]
    A0[:,end,:] .= y_bdy_condi[2]
    A0[:,:,1] .= z_bdy_condi[1]
    A0[:,:,end] .= z_bdy_condi[2]

    B0[1,:,:] .= x_bdy_condi[3]
    B0[end,:,:] .= x_bdy_condi[4]
    B0[:,1,:] .= y_bdy_condi[3]
    B0[:,end,:] .= y_bdy_condi[4]
    B0[:,:,1] .= z_bdy_condi[3]
    B0[:,:,end] .= z_bdy_condi[4]

    C0[1,:,:] .= x_bdy_condi[5]
    C0[end,:,:] .= x_bdy_condi[6]
    C0[:,1,:] .= y_bdy_condi[5]
    C0[:,end,:] .= y_bdy_condi[6]
    C0[:,:,1] .= z_bdy_condi[5]
    C0[:,:,end] .= z_bdy_condi[6]

    # Create the solution containers

    A_container = copy(A0)
    B_container = copy(B0)
    C_container = copy(C0)

    Anew_container = copy(A0)
    Bnew_container = copy(B0)
    Cnew_container = copy(C0)

    total_Em = zeros(i+1)
    total_Em[1] =(1/2)*sum(3*A0.^2)*(δx*δy*δz)    # for keeping track of E_m in the domain
    total_true_Em = copy(total_Em)

    # timestepping loop

    for j in 2:i+1

        display("Currently at timestep j = $(j)")

        for v in 1:iter_steps

            Anew_container = ADI_eqn_A_x!(A_container,Anew_container,B_container,C_container,l_pts,m_pts,n_pts,x_bdy_condi,k1,k2,k3,k4,k5,a6,a7)
            Anew_container = ADI_eqn_A_y!(A_container,Anew_container,B_container,C_container,l_pts,m_pts,n_pts,y_bdy_condi,k1,k2,k3,k4,k5,a6,a7)
            Anew_container = ADI_eqn_A_z!(A_container,Anew_container,B_container,C_container,l_pts,m_pts,n_pts,z_bdy_condi,k1,k2,k3,k4,k5,a6,a7)

            Bnew_container = ADI_eqn_B_x!(A_container,Anew_container,B_container,Bnew_container,C_container,l_pts,m_pts,n_pts,x_bdy_condi,k1,k2,k3,k4,k5,b6,b7,b8,b9)
            Bnew_container = ADI_eqn_B_y!(A_container,Anew_container,B_container,Bnew_container,C_container,l_pts,m_pts,n_pts,y_bdy_condi,k1,k2,k3,k4,k5,b6,b7,b8,b9)
            Bnew_container = ADI_eqn_B_z!(A_container,Anew_container,B_container,Bnew_container,C_container,l_pts,m_pts,n_pts,z_bdy_condi,k1,k2,k3,k4,k5,b6,b7,b8,b9)
            
            Cnew_container = ADI_eqn_C_x!(A_container,Anew_container,B_container,Bnew_container,C_container,Cnew_container,l_pts,m_pts,n_pts,x_bdy_condi,k1,k2,k3,k4,k5,c6,c7,c8)
            Cnew_container = ADI_eqn_C_y!(A_container,Anew_container,B_container,Bnew_container,C_container,Cnew_container,l_pts,m_pts,n_pts,y_bdy_condi,k1,k2,k3,k4,k5,c6,c7,c8)
            Cnew_container = ADI_eqn_C_z!(A_container,Anew_container,B_container,Bnew_container,C_container,Cnew_container,l_pts,m_pts,n_pts,z_bdy_condi,k1,k2,k3,k4,k5,c6,c7,c8)
        
        end

        A_container = copy(Anew_container)
        B_container = copy(Bnew_container)
        C_container = copy(Cnew_container)

        # to keep track of the true solution of A 
        A_true = copy(A0)
        for xco in 1:l_pts
            for yco in 1:m_pts
                for zco in 1:n_pts

                    A_true[xco, yco, zco] = sin((xco-1)*δx)*sin((yco-1)*δy)*sin((zco-1)*δz*(2*pi))*exp(-η*(1+1+(pi/0.5)^2)*(j-1)*δt)
                    
                end
            end
        end

        # so now we have 3 3d arrays for the components of the magnetic field

        # For plotting
        
        if contour=="true"      # B_m is |B|
            if rem(j,10) == 0
                for zcomp in [51] # showing the middle cross sections as contours
                    B_m = zeros(length(x_range), length(y_range))
                    fig_em = Figure()
                    ax_em = Axis(fig_em[1,1], title="z = $(z_range[zcomp]), t = $(δt*(j-1))", xlabel="x", ylabel="y",
                    xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=2))" for value in values ], 
                    ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=2))" for value in yvalues ])
                    for xcomp in 1:length(x_range)
                        for ycomp in 1:length(y_range)
                            B_m[xcomp, ycomp] = sqrt(Anew_container[xcomp,ycomp,zcomp]^2 + Bnew_container[xcomp,ycomp,zcomp]^2 + Cnew_container[xcomp,ycomp,zcomp]^2)
                        end
                    end
                    contour!(ax_em, B_m,levels=50,colormap=:hsv) 
                    Colorbar(fig_em[2,1], vertical=false, limits=(findmin(B_m)[1], findmax(B_m)[1]),colormap=:hsv)
                    display(fig_em)
                end

                for xcomp in [51] 
                    B_m = zeros(length(y_range), length(z_range))
                    fig_em = Figure()
                    ax_em = Axis(fig_em[1,1], title="x = $(round(x_range[xcomp],sigdigits=3)), t = $(δt*(j-1))", xlabel="y", ylabel="z",
                    xtickformat=values -> ["$(round(y_bdy[1]+(value-1)*δy, sigdigits=2))" for value in values ], 
                    ytickformat=yvalues ->["$(round(z_bdy[1]+(value-1)*δz, sigdigits=2))" for value in yvalues ])
                    for zcomp in 1:length(z_range)
                        for ycomp in 1:length(y_range)
                            B_m[ycomp, zcomp] = sqrt(Anew_container[xcomp,ycomp,zcomp]^2 + Bnew_container[xcomp,ycomp,zcomp]^2 + Cnew_container[xcomp,ycomp,zcomp]^2)
                        end
                    end
                    contour!(ax_em, B_m,levels=50,colormap=:hsv) 
                    Colorbar(fig_em[2,1], vertical=false, limits=(findmin(B_m)[1], findmax(B_m)[1]),colormap=:hsv)
                    display(fig_em)
                end

                for ycomp in [51] 
                    B_m = zeros(length(x_range), length(z_range))
                    fig_em = Figure()
                    ax_em = Axis(fig_em[1,1], title="y = $(round(y_range[ycomp],sigdigits=3)), t = $(δt*(j-1))", xlabel="x", ylabel="z",
                    xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=2))" for value in values ], 
                    ytickformat=yvalues ->["$(round(z_bdy[1]+(value-1)*δz, sigdigits=2))" for value in yvalues ])
                    for xcomp in 1:length(x_range)
                        for zcomp in 1:length(z_range)
                            B_m[xcomp, zcomp] = sqrt(Anew_container[xcomp,ycomp,zcomp]^2 + Bnew_container[xcomp,ycomp,zcomp]^2 + Cnew_container[xcomp,ycomp,zcomp]^2)
                        end
                    end
                    contour!(ax_em, B_m,levels=50,colormap=:hsv) 
                    Colorbar(fig_em[2,1], vertical=false, limits=(findmin(B_m)[1], findmax(B_m)[1]),colormap=:hsv)
                    display(fig_em)
                end

            end
        end
        
        if rem(j,10) == 0 # for true soln

            for zcomp in [51] 
                B_m = zeros(length(x_range), length(y_range))
                fig_em = Figure()
                ax_em = Axis(fig_em[1,1], title="z = $(z_range[zcomp]) (true), t = $(δt*(j-1))", xlabel="x", ylabel="y",
                xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=2))" for value in values ], 
                ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=2))" for value in yvalues ])
                for xcomp in 1:length(x_range)
                    for ycomp in 1:length(y_range)
                        B_m[xcomp, ycomp] = sqrt(A_true[xcomp,ycomp,zcomp]^2 + A_true[xcomp,ycomp,zcomp]^2 + A_true[xcomp,ycomp,zcomp]^2)
                    end
                end
                contour!(ax_em, B_m,levels=50,colormap=:hsv) 
                Colorbar(fig_em[2,1], vertical=false, limits=(findmin(B_m)[1], findmax(B_m)[1]),colormap=:hsv)
                display(fig_em)
            end

            for xcomp in [51] 
                B_m = zeros(length(y_range), length(z_range))
                fig_em = Figure()
                ax_em = Axis(fig_em[1,1], title="x = $(round(x_range[xcomp],sigdigits=3)) (true), t = $(δt*(j-1))", xlabel="y", ylabel="z",
                xtickformat=values -> ["$(round(y_bdy[1]+(value-1)*δy, sigdigits=2))" for value in values ], 
                ytickformat=yvalues ->["$(round(z_bdy[1]+(value-1)*δz, sigdigits=2))" for value in yvalues ])
                for zcomp in 1:length(z_range)
                    for ycomp in 1:length(y_range)
                        B_m[ycomp, zcomp] = sqrt(A_true[xcomp,ycomp,zcomp]^2 + A_true[xcomp,ycomp,zcomp]^2 + A_true[xcomp,ycomp,zcomp]^2)
                    end
                end
                contour!(ax_em, B_m,levels=50,colormap=:hsv) 
                Colorbar(fig_em[2,1], vertical=false, limits=(findmin(B_m)[1], findmax(B_m)[1]),colormap=:hsv)
                display(fig_em)
            end

            for ycomp in [51] # Check
                B_m = zeros(length(x_range), length(z_range))
                fig_em = Figure()
                ax_em = Axis(fig_em[1,1], title="y = $(round(y_range[ycomp],sigdigits=3)) (true), t = $(δt*(j-1))", xlabel="x", ylabel="z",
                xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=2))" for value in values ], 
                ytickformat=yvalues ->["$(round(z_bdy[1]+(value-1)*δz, sigdigits=2))" for value in yvalues ])
                for xcomp in 1:length(x_range)
                    for zcomp in 1:length(z_range)
                        B_m[xcomp, zcomp] = sqrt(A_true[xcomp,ycomp,zcomp]^2 + A_true[xcomp,ycomp,zcomp]^2 + A_true[xcomp,ycomp,zcomp]^2)
                    end
                end
                contour!(ax_em, B_m,levels=50,colormap=:hsv) 
                Colorbar(fig_em[2,1], vertical=false, limits=(findmin(B_m)[1], findmax(B_m)[1]),colormap=:hsv)
                display(fig_em)
            end

        end

        # keeping track of the E_m within the domain 
        total_Em[j] = (1/2)*sum(A_container.^2+B_container.^2+C_container.^2)*(δx*δy*δz)
        total_true_Em[j] = (1/2)*sum(3*(A_true.^2))*(δx*δy*δz)

    end

    fig_em = Figure()
    ax_em = Axis(fig_em[1,1], title="E_m within the domain, dotted line is true soln", xlabel="time", ylabel="E_m",
    xtickformat=values -> ["$(round((value-1)*δt, sigdigits=2))" for value in values ])
    lines!(ax_em, total_Em, color=:"blue", linewidth=2)
    lines!(ax_em, total_true_Em, color=:"red", linewidth=3, linestyle=(:dash,:dense))
    display(fig_em)

end
