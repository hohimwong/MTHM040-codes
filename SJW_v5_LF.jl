using CairoMakie
using LinearAlgebra
using Statistics
# code for the JW91 model in spherical coordinate
# ∂A/∂x term treated using leapfrog scheme

function ω(m, δt; amp1 = 0.1, amp2 = 0.1, freq = 1/11)      # non-dimensional function describing the time dependence of omega
Ωᵣ = 1 + amp1*cos(2*pi*freq*(m-1)*δt) + amp2*sin(2*pi*freq*(m-1)*δt)
return Ωᵣ
end

function SJW_v5_LF(param, δq, δt, bdy, m, A0, B0, bdy_condi; a1=0.1, a2=0.1,fq=1/11,EqA_sat_coeff=0,EqB_sat_coeff=0) # bdy_condi is size 4 vector containing bc of A and B respectively
    α0 = copy(param[1]) # this is basically D
    ω0 = copy(param[2])
    r0 = copy(param[3])
    η = copy(param[4])
    λ = copy(param[5])
    n = round(Int, ((bdy[2]-bdy[1])/δq) + 1) # number of grid points
    display(n)
    A0[1] = copy(bdy_condi[1])
    A0[end] = copy(bdy_condi[2])
    B0[1] = copy(bdy_condi[3])
    B0[end] = copy(bdy_condi[4]) # making sure IC satisfies BC

    # now create the coefficient matrices for the equation for A
    γ = η/(r0^2)    # note that γ is just 1 under the non-dimensionalisation
    k = (γ*δt)/(2*δq^2)
    l = (γ*δt)/(2*δq)
    EqA_LHS_matrix = zeros(Float64, n,n)
    q_range = LinRange(bdy[1], bdy[2], n)
    EqA_LHS_matrix[1,1] = 1
    EqA_LHS_matrix[n,n] = 1 
    for i in 2:n-1 
        EqA_LHS_matrix[i,i-1] = k*( (q_range[i])^2 - 1) - l*q_range[i]
        EqA_LHS_matrix[i,i] = 1 + 2*k*(1 - (q_range[i]^2)) + (γ*δt)/(2* (1 - q_range[i]^2) )
        EqA_LHS_matrix[i,i+1] = k*( (q_range[i])^2 - 1) + l*q_range[i]
    end

    EqA_RHS_matrix_1 = zeros(Float64,n,n)
    EqA_RHS_matrix_1[1,1] = 1
    EqA_RHS_matrix_1[n,n] = 1
    for i in 2:n-1
        EqA_RHS_matrix_1[i,i-1] = k*(1 - q_range[i]^2) + l*q_range[i]
        EqA_RHS_matrix_1[i,i] = 1 - 2*k*(1 - q_range[i]^2) - (γ*δt)/(2 * (1 - q_range[i]^2))
        EqA_RHS_matrix_1[i,i+1] = k*(1 - q_range[i]^2) - l*q_range[i]
    end

    EqA_RHS_B_coeff = α0 *δt* q_range 

    # now we create coefficient matrices for the equation for B (note that they are the same as the ones for A)

    EqB_LHS_matrix = zeros(Float64, n,n)
    EqB_LHS_matrix[1,1] = 1
    EqB_LHS_matrix[n,n] = 1
    for i in 2:n-1
        EqB_LHS_matrix[i,i-1] = k*(q_range[i]^2 - 1) - l*q_range[i]
        EqB_LHS_matrix[i,i] = 1 + 2*k*(1 - (q_range[i]^2) ) + (γ*δt)/(2*( 1 - q_range[i]^2 ))
        EqB_LHS_matrix[i,i+1] = l*q_range[i] - k*(1 - q_range[i]^2)
    end

    EqB_RHS_B_matrix = zeros(Float64,n,n)
    EqB_RHS_B_matrix[1,1] = 1
    EqB_RHS_B_matrix[n,n] = 1
    for i in 2:n-1
        EqB_RHS_B_matrix[i,i-1] = k*(1 - q_range[i]^2) + l*q_range[i]
        EqB_RHS_B_matrix[i,i] = 1 - 2*k*(1 - q_range[i]^2) - (γ*δt)/(2*( 1 - q_range[i]^2 ))
        EqB_RHS_B_matrix[i,i+1] = k*(1 - q_range[i]^2) - l*q_range[i]
    end

    EqB_RHS_Anew_matrix = zeros(Float64,n,n) 
    EqB_RHS_Anew_matrix[1,1] = 1
    EqB_RHS_Anew_matrix[n,n] = 1
    EqB_RHS_Aold_matrix = zeros(Float64,n,n)
    EqB_RHS_Aold_matrix[1,1] = 1
    EqB_RHS_Aold_matrix[n,n] = 1

    # now we create the solution containers
    A_list = [zeros(n) for i in 1:m+1] 
    A_list[1] = copy(A0)
    B_list = [zeros(n) for i in 1:m+1]
    B_list[1] = copy(B0)

    # inverting the matrices
    EqA_LHS_matrix_inv = inv(EqA_LHS_matrix)
    EqB_LHS_matrix_inv = inv(EqB_LHS_matrix)

    # timestepping loop
    for i in 2:m+1
        # fill the matrices for the A terms in the B equation
        c1 = r0^2 * ω(i-1,δt, amp1 = a1, amp2 = a2, freq = fq)
        c2 = (c1*δt)/(2)
        c3 = (c1*δt)/(2*2*δq)

        for j in 2:n-1 
            EqB_RHS_Anew_matrix[j,j-1] = c3 * (1 - q_range[j]^2) 
            EqB_RHS_Anew_matrix[j,j] = (c2 * q_range[j]) 
            EqB_RHS_Anew_matrix[j,j+1] = (c3 * ( (q_range[j]^2) -1 ) ) 
        end

        for j in 2:n-1
            EqB_RHS_Aold_matrix[j,j-1] = c3 * (1 - q_range[j]^2) 
            EqB_RHS_Aold_matrix[j,j] = (c2 * q_range[j]) 
            EqB_RHS_Aold_matrix[j,j+1] = (c3 * ( (q_range[j]^2) -1 ) ) # note the two matrices for the A terms in the B equations are the same as well
        end

        A_list[i] = EqA_LHS_matrix_inv*EqA_RHS_matrix_1*A_list[i-1] +
         EqA_LHS_matrix_inv * (EqA_RHS_B_coeff.*B_list[i-1])*(1/(1+EqA_sat_coeff*mean(B_list[i-1].^2))) 
        A_list[i][1] = copy(bdy_condi[1])
        A_list[i][end] = copy(bdy_condi[2]) # applying BC
        B_list[i] = EqB_LHS_matrix_inv * EqB_RHS_B_matrix * B_list[i-1] + 
        (EqB_LHS_matrix_inv * EqB_RHS_Anew_matrix * A_list[i] + EqB_LHS_matrix_inv * EqB_RHS_Aold_matrix * A_list[i-1])*(1/(1+EqB_sat_coeff*mean(B_list[i-1].^2))) - 
         λ*δt*EqB_LHS_matrix_inv*(fill(mean(B_list[i-1].^3), n))
        B_list[i][1] = copy(bdy_condi[3])
        B_list[i][end] = copy(bdy_condi[4]) # applying BC

    end

return A_list, B_list

end
