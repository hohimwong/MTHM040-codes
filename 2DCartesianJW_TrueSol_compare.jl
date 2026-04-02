using CairoMakie
using Statistics
using LinearAlgebra
# code used in checking the diffusion component of the full code
# do note that the "y" here corresponds to "z" in the actual model

# note that in the 2D codes, used m for y, n for x

function ω(iloop, invar; amp=0.5, phase=0.0, freq = 1.7) # function for the time dependence
    output = invar + amp*cos(2*pi*freq*(iloop-1)*δt) + amp*sin(2*pi*freq*(iloop-1)*δt)
    return output
end

# ADI code, the two letters in the end, (A, B) and (x, y) indicates which unknown its for, and which direction it is doing a sweep in
# and again note that the y's here refer to z 

# terms involving ψ will be switched off and set to 0
function adi_eqn_A_y!(container, new_container, B_container, τ, m_pts, n_pts, M1_inv, D, δt, x_range, y_bdy_condi, c1, c2,ψ) # "container" and "new_container" refers to the relevant unknown so here its A
    r = zeros(m_pts)
    B_ms = mean(B_container.^2) 
    for n in 2:n_pts - 1 # z direction
        r[1] = y_bdy_condi[1]
        for m in 2:m_pts-1
            r[m] = c2*container[m,n+1] + c2*container[m,n-1] + c2*new_container[m,n+1] + c2*new_container[m,n-1]+ (c1*c2 - 2*c1 - 2*c2)*container[m,n] + c1*container[m+1,n]+ c1*container[m-1,n] + ψ*δt*D*c1*c2*cos(x_range[n])*B_container[m,n]/(1+τ*B_ms)
        end
        r[end] = y_bdy_condi[2]
        new_container[:,n] = copy(M1_inv*r)
    end
    return new_container
end

function adi_eqn_A_x!(container, new_container, B_container, τ, m_pts, n_pts, M2_inv, D, δt,x_range, x_bdy_condi,c1,c2,ψ)
    s = zeros(n_pts)
    B_ms = mean(B_container.^2)
    for m in 2:m_pts-1 # x direction
        s[1] = x_bdy_condi[1]
        for n in 2:n_pts-1
            s[n] = c2*container[m,n+1] + c2*container[m,n-1] + c1*new_container[m+1,n] + c1*new_container[m-1,n]+ (c1*c2 - 2*c1 - 2*c2)*container[m,n] + c1*container[m+1,n]+ c1*container[m-1,n] + ψ*δt*D*c1*c2*cos(x_range[n])*B_container[m,n]/(1+τ*B_ms)
        end
        s[end] = x_bdy_condi[2]
        new_container[m,:] = copy(M2_inv*s)
    end
    return new_container
end

function adi_eqn_B_y!(container, new_container, A_container, A_new_container , κ, λ, m_pts, n_pts, M1_inv, D, δt, x_range, y_bdy_condi, c1, c2,ψ)
    r = zeros(m_pts)
    B_ms = mean(container.^2)
    B_mc = mean(container.^3)
    for n in 2:n_pts - 1 # z direction
        r[1] = y_bdy_condi[3]
        for m in 2:m_pts-1
            r[m] = c2*container[m,n+1] + c2*container[m,n-1] + c2*new_container[m,n+1] + c2*new_container[m,n-1]+ (c1*c2 - 2*c1 - 2*c2)*container[m,n] + c1*container[m+1,n]+ c1*container[m-1,n] + ψ*( ((δt*sin(x_range[n])*c1*c2)/(2*δx)) * (A_new_container[m,n+1] + A_container[m,n+1])- ((δt*sin(x_range[n])*c1*c2)/(2*δx)) * (A_new_container[m,n] + A_container[m,n]) ) / (1+κ*B_ms) - ψ*c1*c2*λ*δt*B_mc
        end
        r[end] = y_bdy_condi[4]
        new_container[:,n] = copy(M1_inv*r)
    end
    return new_container
end

function adi_eqn_B_x!(container, new_container, A_container, A_new_container, κ, λ, m_pts, n_pts, M2_inv, D, δt,x_range,x_bdy_condi,c1,c2,ψ)
    s = zeros(m_pts)
    B_ms = mean(container.^2)
    B_mc = mean(container.^3)
    for m in 2:m_pts-1 # x direction
        s[1] = x_bdy_condi[3]
        for n in 2:n_pts-1
            s[n] = c2*container[m,n+1] + c2*container[m,n-1] + c1*new_container[m+1,n] + c1*new_container[m-1,n]+ (c1*c2 - 2*c1 - 2*c2)*container[m,n] + c1*container[m+1,n]+ c1*container[m-1,n] + ψ*( ((δt*sin(x_range[n+1])*c1*c2)/(2*δx)) * (A_new_container[m,n+1] + A_container[m,n+1])- ((δt*sin(x_range[n])*c1*c2)/(2*δx)) * (A_new_container[m,n] + A_container[m,n]) ) / (1+κ*B_ms) - ψ*c1*c2*λ*δt*B_mc
        end
        s[end] = x_bdy_condi[4]
        new_container[m,:] = copy(M2_inv*s)
    end
    return new_container
end


function TwoDimCartesianJW_Truesol_comp(param, δx, δy, δt, x_bdy, y_bdy, i, A0, B0, x_bdy_condi, y_bdy_condi; κ=1,τ=1,λ=1,iter_steps=20,source_term_switch=0,update=1)

    c1 = (2*δx^2)/(δt)
    c2 = (2*δy^2)/(δt)

    α0 = copy(param[1])
    ω0 = copy(param[2])
    η = copy(param[3])
    m_pts = round(Int, ((y_bdy[2] - y_bdy[1])/δy) + 1)
    n_pts = round(Int, ((x_bdy[2] - x_bdy[1])/δx) + 1) # grid points in both spatial directions
    x_range = LinRange(x_bdy[1],x_bdy[2],n_pts)

    ψ = 0 # switch for the source terms
    if source_term_switch == 1 # if this = 1 we turn on the source terms
        ψ = 1
    elseif source_term_switch == 0
        ψ = 0
    end

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

    # A and B are matrices now
    
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

    # solution containers for plotting
    A_soln_list = [zeros(m_pts,n_pts) for j in 1:i+1]
    B_soln_list = [zeros(m_pts,n_pts) for j in 1:i+1]
    true_soln_list = [zeros(m_pts,n_pts) for j in 1:i+1]
    A_soln_list[1] = A0
    B_soln_list[1] = B0
    true_soln_list[1] = A0
    err = zeros(n_pts, m_pts)
    max_err_list = zeros(i)
    max_true_soln_list = zeros(i)
    max_numer_soln_list = zeros(i)
    # now we begin time stepping 

    for j in 2:i+1 # timestep loop

        D = α0

        for k in 1:iter_steps # loop for ADI, the auto termination of the loop upon error going under the tolerance level is not implemented here, instead we simply used a large number of iterations that appear to be enough
            Anew_container = adi_eqn_A_y!(A_container, Anew_container, B_container, τ, m_pts, n_pts, M1_inv, D, δt, x_range, y_bdy_condi, c1,c2,ψ)
            Anew_container = adi_eqn_A_x!(A_container, Anew_container, B_container, τ, m_pts, n_pts, M2_inv, D, δt, x_range, x_bdy_condi, c1,c2,ψ)
            Bnew_container = adi_eqn_B_y!(B_container, Bnew_container, A_container, Anew_container, κ, λ, m_pts, n_pts, M1_inv, D, δt, x_range, y_bdy_condi, c1,c2,ψ)
            Bnew_container = adi_eqn_B_x!(B_container, Bnew_container, A_container, Anew_container, κ, λ, m_pts, n_pts, M2_inv, D, δt, x_range, x_bdy_condi, c1,c2,ψ)
        end

        A_container = copy(Anew_container)
        B_container = copy(Bnew_container)
        A_soln_list[j] = copy(A_container)
        B_soln_list[j] = copy(B_container)

        true_soln = zeros(n_pts,m_pts)      # used n_pts for x, m_pts for z here, so need to compare with transpose 
        for xcount in 1:n_pts
            for ycount in 1:m_pts

                true_soln[xcount,ycount] = sin((xcount-1)*δx)*sin(2*(ycount-1)*δy*pi)*exp(-(1+4*pi^2)*(j-1)*δt)

            end
        end
        true_soln_list[j] = copy(true_soln)
        # plot the snapshots of the true and numerical soln at a particular time
        if j == 300
            f_true = Figure()
            ax_true = Axis(f_true[1,1],aspect = 3, title="true solution",xlabel="x",ylabel="z",
            xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=4))" for value in values ], 
            ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=4))" for value in yvalues ])
            ax_num = Axis(f_true[2,1], aspect = 3, title="numerical solution",xlabel="x",ylabel="z",
            xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=4))" for value in values ], 
            ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=4))" for value in yvalues ])
            contour!(ax_true, true_soln, colormap=:hsv,levels=30)
            Colorbar(f_true[1,2],limits=(findmin(true_soln)[1], findmax(true_soln)[1]),colormap=:hsv)
            contour!(ax_num, A_container', colormap=:hsv,levels=30)
            Colorbar(f_true[2,2],limits=(findmin(A_container)[1], findmax(A_container)[1]),colormap=:hsv)
            display(f_true)
        end

        err = abs.(true_soln-A_container) # here we tested A, this (similarly for terms under) can be changed to test for B, where the same results can be observed
        max_err_list[j-1] = findmax(err)[1]
        max_true_soln_list[j-1] = findmax(abs.(true_soln))[1]
        max_numer_soln_list[j-1] = findmax(abs.(A_container))[1]


        # snapshot for B
        if update == 1
            
            if rem(j,1000) == 0 # if j is divisible by 1000 then we make a plot
                fig = Figure()
                axi = Axis(fig[1,1], xlabel="x", ylabel="y", title="Snapshot of B at t = $(δt*(j-1))", 
                xtickformat=values -> ["$(round(x_bdy[1]+(value-1)*δx, sigdigits=4))" for value in values ], 
                ytickformat=yvalues ->["$(round(y_bdy[1]+(value-1)*δy, sigdigits=4))" for value in yvalues ],) 
                contour!(axi, B_soln_list[j]', colormap=:hsv, levels=30)
                Colorbar(fig[1,2], colormap=:hsv, limits=(findmin(B_soln_list[j])[1], findmax(B_soln_list[j])[1]) )
                display(fig)
            end

        end

    end

    true_Em = zeros(length(true_soln_list))
    Em = zeros(length(true_soln_list))

    for s in 1:length(true_soln_list)   # E_m from both the true and numerical soln
        Em[s] = sum((A_soln_list[s].^2)*δx*δy)*(1/2)
        true_Em[s] = sum((true_soln_list[s].^2)*δx*δy)*(1/2)
    end

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="t", ylabel="log(E_m)", title ="Red is true E_m, blue is numerical",
    xtickformat=values -> ["$(round(δt*(value-1), sigdigits=4))" for value in values ])
    lines!(ax, log.(Em), color=:"blue", linewidth=2)
    lines!(ax, log.(true_Em), color=:"red", linewidth=3, linestyle=(:dash,:dense))
    display(fig)

    return A_soln_list, B_soln_list, max_err_list, max_numer_soln_list, max_true_soln_list

end
