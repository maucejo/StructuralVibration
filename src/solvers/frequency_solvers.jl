"""
Structure containing the data feeding the modal solver for calculating an FRF

# Parameters
* ωₙ : Resonance frequencies
* ξₙ : Modal damping ratios
* ϕₑ : Mode shapes at excitation points
* ϕₒ : Mode shapes at observation points

# Note
The mode shapes must be mass-normalized
"""
@with_kw struct ModalFRFProblem
    ωₙ :: Vector{Float64}
    ξₙ :: Vector{Float64}
    ϕₑ :: Matrix{Float64}
    ϕₒ :: Matrix{Float64}

    function ModalFRFProblem(ωₙ, ξₙ, ϕₑ, ϕₒ)
        if !isa(ξₙ, Array)
            ξₙ = fill(ξₙ, length(ωₙ))
        elseif length(ξₙ) != length(ωₙ)
            error("The number of damping ratios must be equal to the number of resonance frequencies")
        end

        new(ωₙ, ξₙ, ϕₑ, ϕₒ)
    end
end

"""
Structure containing the data feeding the direct solver for calculating an FRF

# Parameters
* K : Stiffness matrix
* M : Mass matrix
* C : Damping matrix
* exc_dofs : Degrees of freedom of excitation
* obs_dofs : Degrees of freedom of observation
* freq : Frequencies of interest
"""
@with_kw struct DirectFRFProblem
    K :: Matrix{Float64}
    M :: Matrix{Float64}
    C :: Matrix{Float64}
    exc_dofs :: Vector{Int}
    obs_dofs :: Vector{Int}

    function DirectFRFProblem(K, M, C, exc_dofs, obs_dofs)
        if !isa(exc_dofs, Array)
            exc_dofs = collect(exc_dofs)
        end

        if !isa(obs_dofs, Array)
            obs_dofs = collect(obs_dofs)
        end

        new(K, M, C, exc_dofs, obs_dofs)
    end
end

"""
    solve(m::ModalFRF)
    solve(m::DirectFRF)

Computes the FRF matrix by modal or direct approach

# Parameter
* m : Structure containing the problem data
* freq : Frequencies of interest
* type : Type of FRF to compute (:dis, :vel, :acc)
* ismat : Return the FRF matrix as a 3D array (default = false)

# Output
* FRF : FRF matrix
"""
function solve(m::ModalFRFProblem, freq, type = :dis, ismat = false)
    # Initialisation
    (; ωₙ, ξₙ, ϕₑ, ϕₒ) = m
    Nₑ = size(ϕₑ, 1)
    Nₒ = size(ϕₒ, 1)
    Nf = length(freq)

    FRF = [Matrix{ComplexF64}(undef, Nₒ, Nₑ) for _ in 1:Nf]

    ωf = 2π*freq
    p = Progress(Nf, color = :black, showspeed = true)
    @inbounds for (f, ω) in enumerate(ωf)
        next!(p)
        M = Diagonal(@. 1/(ωₙ^2 - ω^2 + 2im*ξₙ*ωₙ*ω))
        FRF[f] .= ϕₑ*M*ϕₒ'

        if type == :vel
            FRF[f] *= 1im*ω
        elseif type == :acc
            FRF[f] *= -ω^2
        end
    end

    if ismat
        return reshape(reduce(hcat, FRF), Nₒ, Nₑ, :)
    end

    return FRF
end