using CairoMakie
using Statistics
using LinearAlgebra
using SpecialFunctions
include("2D_ang_vel.jl")

# code for the 2D case with our spatial dependence in α and u (ϵ=0)
# do note that the "y" here corresponds to "z" in the actual model, and the code no longer has a switch for source terms !

function u(invar; amp=0.0, phase=0.0, freq = 1.7, l=i, m=m_pts, n=n_pts, x_bdy=x_bdy, y_bdy=y_bdy)       # function taking in ω0 and spits out u with space dependence (3d array)
    
    output = zeros(l+1,m,n)

    spa_dep = TwoD_diff_rot(m,n,x_bdy,y_bdy,1)
    for lloop in 1:l+1
        for mloop in 1:m 
            for nloop in 1:n 
                
                output[lloop,mloop,nloop] = (invar*(1)) * 
                spa_dep[mloop,nloop]
                
            end
        end
    end

    return output
end

function α(invar; amp=0, phase=0.0, freq = 0, l=i, m=m_pts, n=n_pts)        

    output = zeros(l+1,m,n)

    for lloop in 1:l+1
        for mloop in 1:m 
            for nloop in 1:n 

                output[lloop,mloop,nloop] = invar*(1)*
                (cos((nloop-1)*δx)) *
                (1+tanh(4*((mloop-1)*δy-0.4)))
                # the * (1+tanh(4*((mloop-1)*δy-0.4))) bit above can be commented, then we recover the case with only our spatial dependence in u
            end
        end
    end

    return output
    
end

# naming convention for the functions doing the ADI sweeps are the same as before
function adi_eqn_A_y!(container, new_container, B_container, τ, m_pts, n_pts, M1_inv, D, δt, x_range, y_bdy_condi, c1, c2,ψ,α_fn) # "container" and "new_container" refers to the relevant unknown so here its A
    r = zeros(m_pts)
    B_ms = mean(B_container.^2) 
    for n in 2:n_pts - 1 # sweep in z direction
        r[1] = y_bdy_condi[1]
        for m in 2:m_pts-1
            r[m] = c2*container[m,n+1] + c2*container[m,n-1] + c2*new_container[m,n+1] + c2*new_container[m,n-1]+ (c1*c2 - 2*c1 - 2*c2)*container[m,n] + c1*container[m+1,n]+ c1*container[m-1,n] + δt*D*α_fn[m,n]*c1*c2*B_container[m,n]/(1+τ*B_ms)
        end
        r[end] = y_bdy_condi[2]
        new_container[:,n] = copy(M1_inv*r)
    end
    return new_container
end

function adi_eqn_A_x!(container, new_container, B_container, τ, m_pts, n_pts, M2_inv, D, δt,x_range, x_bdy_condi,c1,c2,ψ,α_fn)
    s = zeros(n_pts)
    B_ms = mean(B_container.^2)
    for m in 2:m_pts-1 # x direction
        s[1] = x_bdy_condi[1]
        for n in 2:n_pts-1
            s[n] = c2*container[m,n+1] + c2*container[m,n-1] + c1*new_container[m+1,n] + c1*new_container[m-1,n]+ (c1*c2 - 2*c1 - 2*c2)*container[m,n] + c1*container[m+1,n]+ c1*container[m-1,n] + δt*D*α_fn[m,n]*c1*c2*B_container[m,n]/(1+τ*B_ms)
        end
        s[end] = x_bdy_condi[2]
        new_container[m,:] = copy(M2_inv*s)
    end
    return new_container
end

function adi_eqn_B_y!(container, new_container, A_container, A_new_container , κ, λ, m_pts, n_pts, M1_inv, D, δt, x_range, y_bdy_condi, c1, c2,ψ,α_fn,u_fn,α_fn_new,u_fn_new,δx,δy,ϵ,τ)
    r = zeros(m_pts)
    B_ms = mean(container.^2)
    B_mc = mean(container.^3)
    for n in 2:n_pts - 1 # z direction
        r[1] = y_bdy_condi[3]
        for m in 2:m_pts-1
            r[m] = c2*container[m,n+1] + c2*container[m,n-1] + c2*new_container[m,n+1] + c2*new_container[m,n-1] +
            (c1*c2 - 2*c1 - 2*c2)*container[m,n] + c1*container[m+1,n]+ c1*container[m-1,n] +
            ( (δt/2)*(1/(1+κ*B_ms))*( ((u_fn_new[m+1,n]-u_fn_new[m-1,n])/(2*δy))*((A_new_container[m,n+1]-A_new_container[m,n-1])/(2*δx)) + ((u_fn[m+1,n]-u_fn[m-1,n])/(2*δy))*((A_container[m,n+1]-A_container[m,n-1])/(2*δx))) -
            (ψ*δt/2)*(1/(1+κ*B_ms))*( ((u_fn_new[m,n+1]-u_fn_new[m,n-1])/(2*δx))*((A_new_container[m+1,n]-A_new_container[m-1,n])/(2*δy)) + ((u_fn[m,n+1]-u_fn[m,n-1])/(2*δx))*((A_container[m+1,n]-A_container[m-1,n])/(2*δy))) -
            (ϵ*δt/2)*(1/(1+τ*B_ms))*( ((α_fn_new[m+1,n]-α_fn_new[m-1,n])/(2*δy))*((A_new_container[m+1,n]-A_new_container[m-1,n])/(2*δy)) + ((α_fn_new[m,n+1]-α_fn_new[m,n-1])/(2*δx))*((A_new_container[m,n+1]-A_new_container[m,n-1])/(2*δx)) +
            ((α_fn[m+1,n]-α_fn[m-1,n])/(2*δy))*((A_container[m+1,n]-A_container[m-1,n])/(2*δy)) + ((α_fn[m,n+1]-α_fn[m,n-1])/(2*δx))*((A_container[m,n+1]-A_container[m,n-1])/(2*δx)) +
            (α_fn_new[m,n]*(A_new_container[m,n+1]-2*A_new_container[m,n]+A_new_container[m,n-1])/(δx^2)) +
            (α_fn_new[m,n]*(A_new_container[m+1,n]-2*A_new_container[m,n]+A_new_container[m-1,n])/(δy^2)) +
            (α_fn[m,n]*(A_container[m,n+1]-2*A_container[m,n]+A_container[m,n-1])/(δx^2)) +
            (α_fn[m,n]*(A_container[m+1,n]-2*A_container[m,n]+A_container[m-1,n])/(δy^2))) ) * (c1*c2) - c1*c2*λ*δt*B_mc
        end
        r[end] = y_bdy_condi[4]
        new_container[:,n] = copy(M1_inv*r)
    end
    return new_container
end

function adi_eqn_B_x!(container, new_container, A_container, A_new_container, κ, λ, m_pts, n_pts, M2_inv, D, δt,x_range,x_bdy_condi,c1,c2,ψ,α_fn,u_fn,α_fn_new,u_fn_new,δx,δy,ϵ,τ)
    s = zeros(m_pts)
    B_ms = mean(container.^2)
    B_mc = mean(container.^3)
    for m in 2:m_pts-1 # x direction
        s[1] = x_bdy_condi[3]
        for n in 2:n_pts-1
            s[n] = c2*container[m,n+1] + c2*container[m,n-1] + c1*new_container[m+1,n] + c1*new_container[m-1,n] +
            (c1*c2 - 2*c1 - 2*c2)*container[m,n] + c1*container[m+1,n]+ c1*container[m-1,n] + 
            ( (δt/2)*(1/(1+κ*B_ms))*( ((u_fn_new[m+1,n]-u_fn_new[m-1,n])/(2*δy))*((A_new_container[m,n+1]-A_new_container[m,n-1])/(2*δx)) + ((u_fn[m+1,n]-u_fn[m-1,n])/(2*δy))*((A_container[m,n+1]-A_container[m,n-1])/(2*δx))) -
            (ψ*δt/2)*(1/(1+κ*B_ms))*( ((u_fn_new[m,n+1]-u_fn_new[m,n-1])/(2*δx))*((A_new_container[m+1,n]-A_new_container[m-1,n])/(2*δy)) + ((u_fn[m,n+1]-u_fn[m,n-1])/(2*δx))*((A_container[m+1,n]-A_container[m-1,n])/(2*δy))) -
            (ϵ*δt/2)*(1/(1+τ*B_ms))*( ((α_fn_new[m+1,n]-α_fn_new[m-1,n])/(2*δy))*((A_new_container[m+1,n]-A_new_container[m-1,n])/(2*δy)) + ((α_fn_new[m,n+1]-α_fn_new[m,n-1])/(2*δx))*((A_new_container[m,n+1]-A_new_container[m,n-1])/(2*δx)) +
            ((α_fn[m+1,n]-α_fn[m-1,n])/(2*δy))*((A_container[m+1,n]-A_container[m-1,n])/(2*δy)) + ((α_fn[m,n+1]-α_fn[m,n-1])/(2*δx))*((A_container[m,n+1]-A_container[m,n-1])/(2*δx)) +
            (α_fn_new[m,n]*(A_new_container[m,n+1]-2*A_new_container[m,n]+A_new_container[m,n-1])/(δx^2)) +
            (α_fn_new[m,n]*(A_new_container[m+1,n]-2*A_new_container[m,n]+A_new_container[m-1,n])/(δy^2)) +
            (α_fn[m,n]*(A_container[m,n+1]-2*A_container[m,n]+A_container[m,n-1])/(δx^2)) +
            (α_fn[m,n]*(A_container[m+1,n]-2*A_container[m,n]+A_container[m-1,n])/(δy^2))) ) * (c1*c2) - c1*c2*λ*δt*B_mc
        end
        s[end] = x_bdy_condi[4]
        new_container[m,:] = copy(M2_inv*s)
    end
    return new_container
end


function TwoDJW_alp2ome_v3(param, δx, δy, δt, x_bdy, y_bdy, i, A0, B0, x_bdy_condi, y_bdy_condi; κ=1,τ=1,λ=1,ux=1,update=1,ϵ=0,tol=10e-2)

    c1 = (2*δx^2)/(δt)
    c2 = (2*δy^2)/(δt)

    α0 = copy(param[1])
    ω0 = copy(param[2])
    η = copy(param[3])
    L = copy(param[4])
    m_pts = round(Int, ((y_bdy[2] - y_bdy[1])/δy) + 1)
    n_pts = round(Int, ((x_bdy[2] - x_bdy[1])/δx) + 1)
    x_range = LinRange(x_bdy[1],x_bdy[2],n_pts)

    ψ = 0 
    if ux == 1 # if this = 1 we turn on the x derivative of u
        ψ = 1
    elseif ux == 0
        ψ = 0
    end

    α_container = α(1)
    u_container = u(ω0)
    

    # create coefficient matrix
    M1 = zeros(m_pts,n_pts)
    M1[1,1] = 1
    M1[end,end] = 1
    for j in 2:m_pts-1
        M1[j,j-1] = -c1
        M1[j,j] = c1*c2 + 2*c1 + 2*c2
        M1[j,j+1] = -c1
    end
    M1_inv = inv(M1)

    M2 = zeros(m_pts,n_pts)
    M2[1,1] = 1
    M2[end,end] = 1
    for j in 2:m_pts-1
        M2[j,j-1] = -c2
        M2[j,j] = c1*c2 + 2*c1 + 2*c2
        M2[j,j+1] = -c2
    end
    M2_inv = inv(M2)

    # A and B are matrices
    
    # making sure the BC is applied to the initial state
    A0[1,:] .= y_bdy_condi[1]
    A0[end,:] .= y_bdy_condi[2]
    A0[:,1] .= x_bdy_condi[1]
    A0[:,end] .= x_bdy_condi[2]
    
    B0[1,:] .= y_bdy_condi[3]
    B0[end,:] .= y_bdy_condi[4]
    B0[:,1] .= x_bdy_condi[3]
    B0[:,end] .= x_bdy_condi[4]

    # setting up the initial state 

    A_container = copy(A0)
    B_container = copy(B0)

    Anew_container = copy(A_container)
    Bnew_container = copy(B_container)

    A_prev = copy(A0)
    B_prev = copy(B0)

    # solution containers for plotting
    A_soln_list = [zeros(m_pts,n_pts) for j in 1:i+1]
    B_soln_list = [zeros(m_pts,n_pts) for j in 1:i+1]
    A_soln_list[1] = A0
    B_soln_list[1] = B0
    D = (α0*ω0*(L^3))/(η^2)
    @show(D)
    # now we begin time stepping 

    for j in 2:i+1 # timestep loop

        if rem(j,10)==0
            display("Currently at j = $(j)")
        end
        
        err_A = 100.0     # arbitrary number that is larger than the error tolerance
        err_B = 100.0

        u_fn = u_container[j-1,:,:]     # getting u and α at the given time
        α_fn = α_container[j-1,:,:]

        u_fn_new = u_container[j,:,:]
        α_fn_new = α_container[j,:,:]

        while >(err_A,tol) || >(err_B,tol) # loop for adi, stops when converged sufficiently

            Anew_container = adi_eqn_A_y!(A_container, Anew_container, B_container, τ, m_pts, n_pts, M1_inv, D, δt, x_range, y_bdy_condi, c1,c2,ψ,α_fn)
            Anew_container = adi_eqn_A_x!(A_container, Anew_container, B_container, τ, m_pts, n_pts, M2_inv, D, δt, x_range, x_bdy_condi, c1,c2,ψ,α_fn)
            Bnew_container = adi_eqn_B_y!(B_container, Bnew_container, A_container, Anew_container, κ, λ, m_pts, n_pts, M1_inv, D, δt, x_range, y_bdy_condi, c1,c2,ψ,α_fn,u_fn,α_fn_new,u_fn_new,δx,δy,ϵ,τ)
            Bnew_container = adi_eqn_B_x!(B_container, Bnew_container, A_container, Anew_container, κ, λ, m_pts, n_pts, M2_inv, D, δt, x_range, x_bdy_condi, c1,c2,ψ,α_fn,u_fn,α_fn_new,u_fn_new,δx,δy,ϵ,τ)
            
            err_A = findmax(abs.(Anew_container-A_prev))[1]
            err_B = findmax(abs.(Bnew_container-B_prev))[1]
            A_prev = copy(Anew_container)
            B_prev = copy(Bnew_container)
        
        end

        A_container = copy(Anew_container)
        B_container = copy(Bnew_container)
        A_soln_list[j] = copy(A_container)
        B_soln_list[j] = copy(B_container)

        if update == 1
            
            if rem(j,5) == 0 # if j is divisible by the second input then we make a plot
                fig = Figure()
                axi1 = Axis(fig[1,1], xlabel="x", ylabel="z", title="Snapshot of B at t = $(round(δt*(j-1),digits=3))", 
                xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=4))" for value in values ], 
                ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=4))" for value in yvalues ], aspect=2.5) 
                contour!(axi1, B_soln_list[j]', colormap=:hsv, levels=30)
                Colorbar(fig[1,2], colormap=:hsv, limits=(findmin(B_soln_list[j])[1], findmax(B_soln_list[j])[1]) )

                axi2 = Axis(fig[2,1], xlabel="x", ylabel="z", title="Snapshot of A at t = $(round(δt*(j-1),digits=3))", 
                xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=4))" for value in values ], 
                ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=4))" for value in yvalues ], aspect=2.5) 
                contour!(axi2, A_soln_list[j]', colormap=:viridis, levels=30)
                Colorbar(fig[2,2], colormap=:viridis, limits=(findmin(A_soln_list[j])[1], findmax(A_soln_list[j])[1]) )
                display(fig)
            end

        end

    end
    return A_soln_list, B_soln_list

end