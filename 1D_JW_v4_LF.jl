using CairoMakie
using LinearAlgebra
using Statistics

# nondimensionalised 1D Cartesian JW code
# uses leapfrog for the A derivitive in the ∂B/∂t equation

# non-dimensional function for time dep.
function ω(m, δt; amp1 = 0.0, amp2 = 0.0, freq = 0.2) # adjust the amplitude/frequency of the time dependence here by changing amp1 amp2 (amplitude) or freq (frequency)
Ωᵣ = 1 + amp1*cos(2*pi*freq*(m-1)*δt) + amp2*sin(2*pi*freq*(m-1)*δt)
return Ωᵣ
end

function OneD_JW_v4_LF(param, δx, δt, bdy, m, A0, B0) # bdy is boundary, size 2 vector 
    # param = [α0, ω0, r0, η, κ, τ, λ] A0 B0 initial values
    α0 = param[1] # α0 is effectively the dynamo number in our non-dimensionalisation
    ω0 = param[2]
    r0 = param[3]
    η = param[4]
    κ = param[5]
    τ = param[6]
    λ = param[7]
    k = δt/(2*(δx)^2)
    n = round(Int, ((bdy[2]-bdy[1])/δx) + 1) # number of grid points
    A0[1] = 0
    A0[end] = 0
    B0[1] = 0
    B0[end] = 0 # making sure IC satisfies BC
    x_range = LinRange(bdy[1],bdy[2],n) # the grid

    # creating the matrices containing the coefficients
    T1 = zeros(Float64, n, n)
    T1[1,1] = 1
    T1[n,n] = 1
    for i in 2:n-1 # filling in the elements of the matrix T₁ (LHS matrix)
        T1[i,i-1] = -k
        T1[i,i] = 1+2*k
        T1[i, i+1] = -k
    end
    T2 = zeros(Float64,n,n) # RHS matrix
    T2[1,1] = 1 
    T2[n,n] = 1 
    for i in 2:n-1
        T2[i,i-1] = k
        T2[i,i] = 1 - 2*k
        T2[i,i+1] = k
    end
    M1 = zeros(Float64, n, n) # same as T1
    M1[1,1] = 1
    M1[n,n] = 1
    for i in 2:n-1
        M1[i,i-1] = -k 
        M1[i,i] = 1+2*k 
        M1[i,i+1] = -k
    end
    M2 = zeros(Float64,n,n) # same as T2
    M2[1,1] = 1 
    M2[n,n] = 1 
    for i in 2:n-1
        M2[i,i-1] = k 
        M2[i,i] = 1-2*k 
        M2[i,i+1] = k 
    end
    
    A_list = [zeros(n) for i in 1:m+1] # creating solution containers
    A_list[1] = copy(A0)
    B_list = [zeros(n) for i in 1:m+1]
    B_list[1] = copy(B0)
    T1_inverse = copy(inv(T1))
    RHSmatrix1 = T1_inverse*T2 # for the equation for A 
    M1_inverse = copy(inv(M1))
    RHSmatrix2 = M1_inverse*M2 # for the equation for B
    EqA_inhomog_coeff = zeros(n)
    EqB_inhomog_matrix = zeros(Float64,n,n) # for the inhomogeous term in equation for B 
    EqB_inhomog_matrix[1,1] = 1
    EqB_inhomog_matrix[n,n] = 1
    for i in 2:n-1
        EqB_inhomog_matrix[i,i-1] = -1*sin(x_range[i])
        EqB_inhomog_matrix[i,i+1] = 1*sin(x_range[i])
    end

    D = (α0*ω0*(r0)^3)/(η^2)    # every parameter except α is set to 1
    EqA_inhomog_coeff = [D*δt*cos(x_range[j]) for j in 1:n]

    for i in 2:m+1      # timestepping loop
        EqB_inhomog_coeff = (δt/(4*δx))*ω(i-1,δt)
        A_list[i] = RHSmatrix1*(A_list[i-1]) + T1_inverse*( ( 1/(1+τ*mean(B_list[i-1].^2) ) ) * EqA_inhomog_coeff.*B_list[i-1])
        A_list[i][1] = 0 # applying BC
        A_list[i][end] = 0  
        B_list[i] = RHSmatrix2*(B_list[i-1]) + (M1_inverse*EqB_inhomog_matrix)*((1/(1+κ*mean(B_list[i-1].^2)))*EqB_inhomog_coeff*A_list[i-1]) +
         (M1_inverse*EqB_inhomog_matrix)*((1/(1+κ*mean(B_list[i-1].^2)))*EqB_inhomog_coeff*A_list[i]) -
          λ*δt*(M1_inverse*fill(( (mean(B_list[i-1].^3))),n) ) 
        # note we have an extra term (scalar) encapsulating the coefficients/time dep.; this is the effectively the same as what's written in the report
        B_list[i][1] = 0 # applying BC
        B_list[i][end] = 0  
    end

return A_list, B_list
end

