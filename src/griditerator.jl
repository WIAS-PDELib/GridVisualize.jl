struct LinearSimplices{D, Tc, Ti, Tf} <: LinearSimplexIterator{D}
    coord::Matrix{Tc}
    cellnodes::Matrix{Ti}
    values::Vector{Tf}
    gridscale::Tc
    range::StepRange{Int, Int}
end

function LinearSimplices(coord::Matrix{Tc}, cn::Matrix{Ti}, f::Vector{Tf}; gridscale = 1.0, nthreads = Threads.nthreads()) where {Tc, Ti, Tf}
    ncells = size(cn, 2)
    dim = size(coord, 1)
    return map(enumerate(chunks(1:ncells; n = nthreads))) do c
        LinearSimplices{dim, Tc, Ti, Tf}(coord, cn, f, gridscale, last(c))
    end
end

function LinearSimplices(g::ExtendableGrid, f::Vector; nthreads = Threads.nthreads(), gridscale = 1.0)
    return LinearSimplices(g[Coordinates], g[CellNodes], f; nthreads, gridscale)
end

function Base.iterate(linear_simplices::LinearSimplices{D}, args...) where {D}
    (; coord, cellnodes, values, gridscale, range) = linear_simplices
    iter = iterate(range, args...)
    isnothing(iter) && return nothing
    (icell, state) = iter
    @views s = LinearSimplex(
        Val{D},
        coord[:, cellnodes[:, icell]],
        values[cellnodes[:, icell]],
        gridscale
    )
    return (s, state)
end
