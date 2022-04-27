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
    return copysign(sqrt(g2test(m)), m[1,1] - expected)
end

md"""
Computes the root of the G^2 test statistic for a 2x2 table presented as
four separate entries that are arranged into a contingency table this way
```
k11 k12
k21 k22
```
"""
signedG2(k11, k12, k21, k22) = signedG2([k11 k12; k21 k22])

md"""
Computes signed G statistics for the row-wise cooccurence of 
non-zero elements in a matrix.

This matrix is arranged with items in columns and users/events/windows 
in rows. 
"""
indicators(A::Matrix; kw...) = indicators(sparse(A), kw...)

function indicators(A; itemcut=200)
    # clip all non-zero elements of A to 1
    rows, items = size(A)
    for (i, j, v) in zip(findnz(A)...)
        A[i, j] = 1
    end

    if itemcut ≤ 0
        itemcut = rows
    end

    # force each column to have a limited number of non-zeros
    for i in 1:items
        jx = findnz(A[:, i])[1]
        excess = length(jx) - itemcut
        if excess > 0
            drops = shuffle!(jx)[1:excess]
            A[drops, i] .= 0
        end
    end
    dropzeros!(A)
    
    # and compute item (column) totals
    itemCounts = [sum(A[:,i]) for i in 1:items]
    total = sum(itemCounts)
    
    # compute cooccurrence
    cooc = A' * A

    # and scores
    nonzeros = findnz(cooc)
    ix = Vector{Int}()
    jx = Vector{Int}()
    rx = Vector{Float64}()
    for (i, j, v) in zip(nonzeros...)
        if i < j
            k11 = cooc[i, j]
            k12 = itemCounts[i] - k11
            k21 = itemCounts[j] - k11
            k22 = rows - k11 - k12 - k21
            score = signedG2(k11, k12, k21, k22)
            push!(ix, i)
            push!(jx, j)
            push!(rx, score)

            push!(ix, j)
            push!(jx, i)
            push!(rx, score)
        end
    end
    sparse(ix, jx, rx)
end


function denormEntropy(v)
    N = sum(v)
    sum(-v .* log.(v/N + (v .== 0)))
end

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
