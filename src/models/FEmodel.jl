"""
    Mesh(s, xmin, L, Nelt)

Construct a mesh for a beam with Nelt elements, length L and starting at xmin.

# Inputs
- s: Structure containing the data related to the 1D system
- xmin: starting position of the beam
- L: length of the beam
- Nelt: number of elements
- bc: Boundary conditions type
    * :CC : Clamped - Clamped
    * :FF : Free - Free
    * :CF : Clamped - Free
    * :SS : Simply Supported - Simply Supported (specific to beam)
    * :CS : Clamped - Simply Supported (specific to beam)
    * :SF : Simply Supported - Free (specific to beam)


# Outputs
- mesh: mesh of the beam

# Example
```julia-repl
julia> mesh = Mesh(0., 1., 10)
```
"""
@with_kw struct Mesh
    xmin :: Float64
    L :: Float64
    Nelt :: Int64
    Nnodes :: Int64
    Nodes :: Matrix{Float64}
    Elt :: Matrix{Int64}
    Ndof_per_node :: Int64
    elem_size :: Float64
    constrained_dofs :: Vector{Int64}
    free_dofs :: Vector{Int64}

    function Mesh(s :: BeamBarRod, xmin, L, Nelt, bc = :CC)
        Nnodes = Nelt + 1
        Nodes = zeros(Nnodes, 2)
        Elt = zeros(Nelt, 3)
        elem_size = L/Nelt

        for i = 1:Nnodes
            Nodes[i, 1] = i
            Nodes[i, 2] = xmin + (i - 1)*elem_size
        end

        for i = 1:Nelt
            Elt[i, 1] = i
            Elt[i, 2] = i
            Elt[i, 3] = i + 1
        end

        if isa(s, Beam)
            Ndof_per_node = 2
            dofs = collect(1:2Nnodes)
            if bc == :SS
                constrained_dofs = [1, 2Nnodes - 1]
            elseif bc == :CC
                constrained_dofs = [1, 2, 2Nnodes - 1, 2Nnodes]
            elseif bc == :CS
                constrained_dofs = [1, 2, 2Nnodes - 1]
            elseif bc == :CF
                constrained_dofs = [1, 2]
            elseif bc == :SF
                constrained_dofs = [1]
            elseif bc == :FF
                constrained_dofs = Int[]
            else
                error("Boundary conditions not implemented")
            end
        elseif isa(s, BarRod)
            Ndof_per_node = 1
            dofs = collect(1:Nnodes)
            if bc == :CC
                constrained_dofs = [1, Nnodes]
            elseif bc == :CF
                constrained_dofs = [1]
            elseif bc == :FF
                constrained_dofs = Int[]
            else
                error("Boundary conditions not implemented")
            end
        end
        free_dofs = setdiff(dofs, constrained_dofs)

        new(xmin, L, Nelt, Nnodes, Nodes, Elt, Ndof_per_node, elem_size, constrained_dofs, free_dofs)
    end
end

"""
    assembly(s, mesh)

Compute the global stiffness and mass matrices for a beam with a given mesh.

# Inputs
- s: Structure containing the data related to the 1D system
- mesh: Mesh

# Outputs
- `K`: global stiffness matrix
- `M`: global mass matrix

# Example
```julia-repl
julia> b = Beam(1., 3e-2, 1e-2, 2.1e11, 7850.)
julia> mesh = BeamMesh(0., 1., 10)
julia> K, M = assembly(b, mesh)
```
"""
function assembly(s :: BeamBarRod, mesh :: Mesh)
    # Compute elemental matrices
    kₑ, mₑ = element_matrix(s, mesh.elem_size)

    (; Nnodes, Nelt, Ndof_per_node) = mesh

    # Assemble global matrices
    Nddl = Nnodes*Ndof_per_node

    K = zeros(Nddl, Nddl)
    M = zeros(Nddl, Nddl)
    ind = zeros(eltype(Nnodes), Ndof_per_node^2)
    @inbounds @views for i = 1:Nelt
        ind .= (Ndof_per_node*mesh.Elt[i, 2:end]' .+ repeat((0:Ndof_per_node-1), 1, Ndof_per_node) .- Ndof_per_node .+ 1)[:];

        K[ind, ind] += kₑ
        M[ind, ind] += mₑ
    end

    return K, M
end

"""
    element_matrix(s::Beam, h)
    element_matrix(s::BarRod, h)

Compute the elemental stiffness and mass matrices for a beam with a element size `h`.

# Inputs
- s: Structure containing the data related to the 1D system
- h: element size

# Outputs
- `kₑ`: elemental stiffness matrix
- `mₑ`: elemental mass matrix

# Example
```julia-repl
julia> b = Beam(1., 3e-2, 1e-2, 2.1e11, 7850.)
julia> kₑ, mₑ = beam_elem_elt(b, 1e-2)
```
"""
function element_matrix(beam :: Beam, h)
    # Constants
    kc = beam.D/h^3
    mc = beam.m*h/420.

    # Elemental stiffness matrix
    kₑ = kc.*[12. 6h -12. 6h;
              6h 4h^2 -6h 2h^2;
              -12. -6h 12. -6h;
              6h 2h^2 -6h 4h^2]

    # Elemental mass matrix
    mₑ = mc.*[156. 22h 54. -13h;
              22h 4h^2 13h -3h^2;
              54. 13h 156. -22h;
              -13h -3h^2 -22h 4h^2]

    return kₑ, mₑ
end

function element_matrix(barod :: BarRod, h)
    # Constants
    kc = barod.D/h
    mc = barod.m*h/6.

    # Elemental stiffness matrix
    kₑ = kc.*[1. -1.;
              -1. 1.]

    # Elemental mass matrix
    mₑ = mc.*[2. 1.;
              1. 2.]

    return kₑ, mₑ
end

"""
    dofs_selection(mesh, X)

Select the dofs corresponding to the closest nodes to the positions X.

# Inputs
- mesh: Structure mesh
- `X`: positions

# Outputs
- `dofs`: dofs corresponding to the closest nodes
- `S`: selection matrix

# Example
```julia-repl
julia> mesh = BeamMesh(0., 1., 10)
julia> dofs, S = dofs_selection(mesh, [0.1, 0.2])
```
"""
function dofs_selection(mesh :: Mesh, X)
    N = length(X)
    dofs = zeros(Int, N)
    S = zeros(N, length(mesh.free_dofs))
    @inbounds for (i, Xi) in enumerate(X)
        d = @. abs(mesh.Nodes[:, 2] - Xi)
        if mesh.Ndof_per_node == 2
            dofs[i] = 2argmin(d) - 1
        else
            dofs[i] = argmin(d)
        end

        pos = findall(mesh.free_dofs .== dofs[i])
        if length(pos) != 0
            S[i, pos[1]] = 1.
        end
    end

    return dofs, S
end