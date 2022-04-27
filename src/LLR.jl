module LLR

using Markdown
using SparseArrays
using Random

md"""
Computes the G^2 test statistic for a contingency table in a 2x2 matrix

This test is a \chi^2 test that avoids the assumption of normal distribution
inherent in Pearson's \chi^2 test. This is very helpful for sparse problems.

See http://tdunning.blogspot.com/2008/03/surprise-and-coincidence.html and https://www.researchgate.net/publication/2477641_Accurate_Methods_for_the_Statistics_of_Surprise_and_Coincidence for more information.
"""
function g2test(m::Matrix)::Float64
    2 * (-denormEntropy(m) + denormEntropy(sum(m, dims=1)) + denormEntropy(sum(m, dims=2)))
end

md"""
Computes the root of the G^2 test statistic for a 2x2 table.

For uncorrelated rows and columns and large enough counts, this should
be normally distributed.
"""
function signedG2(m::Matrix)::Float64
    if size(m) != (2, 2)
        raise(ArgumentError("Must have 2x2 matrix"))
    end
    total = sum(m)
    expected = sum(m[:,1]) / total * sum(m[1,:])
    return copysign(sqrt(max(0, g2test(m))), m[1,1] - expected)
end

md"""
Computes the root of the G^2 test statistic for a 2x2 table presented as
four separate entries that are arranged into a contingency table this way
```
k11 k12
k21 k22
```
This version, however, is unwound to avoid allocating small matrices
"""
function signedG2(k11, k12, k21, k22)
    full = denormEntropy(k11, k12, k21, k22)
    row = denormEntropy(k11 + k12, k21 + k22)
    col = denormEntropy(k11 + k21, k12 + k22)
    expected = (k11 + k21) / (k11 + k12 + k21 + k22) * (k11 + k12)
    copysign(sqrt(max(0, 2 * (full - row - col))), k11 - expected)
end

md"""
Computes signed G statistics for the row-wise cooccurence of 
non-zero elements in a matrix.

This matrix is arranged with items in columns and users/events/windows 
in rows. 
"""
indicators(A::Matrix; kw...) = indicators(sparse(A), kw...)

function indicators(A; itemcut=0, rowcut=200)
    # clip all non-zero elements of A to 1
    rows, items = size(A)
    for (i, j, v) in zip(findnz(A)...)
        A[i, j] = 1
    end

    if itemcut ≤ 0
        itemcut = rows
    end

    # force each row to have a limited number of non-zeros
    # this has minimal impact on the recommendations and
    # makes the computation very fast 
    if rowcut > 0
        for i in 1:rows
            jx = findnz(A[i, :])[1]
            excess = length(jx) - rowcut
            if excess > 0
                drops = shuffle!(jx)[1:excess]
                A[i, drops] .= 0
            end
        end
    end

    # likewise for each column 
    # (this is disabled by default because it has little impact)
    if itemcut > 0
        for j in 1:items
            ix = findnz(A[:, j])[1]
            excess = length(ix) - itemcut
            if excess > 0
                drops = shuffle!(ix)[1:excess]
                A[drops, j] .= 0
            end
        end
    end
    dropzeros!(A)
    @info "cut complete"
    
    # and compute item (column) totals
    itemCounts = [sum(A[:,i]) for i in 1:items]
    total = sum(itemCounts)
    
    # compute cooccurrence
    cooc = A' * A
    @info "cooc done"

    # and scores
    nonzeros = findnz(cooc)
    ix = Vector{Int}()
    jx = Vector{Int}()
    k11 = Vector{Float64}()
    k1x = Vector{Float64}()
    kx1 = Vector{Float64}()
    for (i, j, v) in zip(nonzeros...)
        if i < j
            push!(k11, cooc[i, j])
            push!(k1x, itemCounts[i])
            push!(kx1, itemCounts[j])
            push!(ix, i)
            push!(jx, j)

            push!(k11, cooc[i, j])
            push!(k1x, itemCounts[i])
            push!(kx1, itemCounts[j])
            push!(ix, j)
            push!(jx, i)
        end
    end
    k12 = k1x .- k11
    k21 = kx1 .- k11
    k22 = rows .- (k11 .+ k12 .+ k21)
    @info "setup complete"
    sparse(ix, jx, signedG2.(k11, k12, k21, k22))
end


function denormEntropy(v)
    N = sum(v)
    sum(-v .* log.(v/N + (v .== 0)))
end

md"This version doesn't allocate and is easier to inline"
function denormEntropy(k1, k2)
    N = k1 + k2
    kLogP(k1, k1 / N) + kLogP(k2, k2 / N)
end

md"This version doesn't allocate and is easier to inline"
function denormEntropy(k1, k2, k3, k4)
    N = k1 + k2 + k3 + k4
    kLogP(k1, k1 / N) + kLogP(k2, k2 / N) + 
        kLogP(k3, k3 / N) + kLogP(k4, k4 / N)
end

kLogP(k, p) = k * log(p + (p == 0))

"""
Compares the frequencies of corresponding items in two dictionaries using LLR
"""
function compare(a, b)
    aTotal = sum(values(a))
    bTotal = sum(values(b))
    let f = key -> 
        let va = get(a, key, 0), vb = get(b, key, 0) 
            return (key => signedG2(va, aTotal - va, vb, bTotal - vb))
        end
        return Dict(f(key) for key in union(keys(a), keys(b)))
    end
end


end # module
