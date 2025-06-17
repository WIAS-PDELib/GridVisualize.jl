using ExtendableGrids
using GridVisualizeTools
using GridVisualize
using Test


function testloops(dim)
    X = 0:0.1:10
    if dim == 1
        g = simplexgrid(X)
    elseif dim == 2
        g = simplexgrid(X, X)
    else
        g = simplexgrid(X, X, X)
    end
    f = ones(num_nodes(g))
    ls = LinearSimplices(g, f; nthreads = 3)
    testloop(ls) # for compilation
    nalloc = @allocated sum_f = testloop(ls)


    @test nalloc < 256 # allow for some allocations
    return @test sum_f == (dim + 1) * num_cells(g)
end

@testset "iterator testloops" begin
    testloops(1)
    testloops(2)
    testloops(3)
end
