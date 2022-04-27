# LLR
This package is a Julia native implementation of the G^2 test in Julia
along with related routines for comparing frequencies and cooccurrences.

For more information on the original G^2 usage in corpus linguistics,
see the 1993 [Surprise and Coincidence](https://aclanthology.org/J93-1003.pdf) paper.

For more information on how to use this test for recommendations, see
this paper on [Efficient Incremental Cooccurrence Analysis](https://ssc.io/pdf/p3-schelter.pdf).

# Usage for Recommendations
If you have a matrix `A` that contains a history of interactions between users and items where
rows represent users and columns represent items, you can reduce that matrix to one that records
which items are interestingly cooccurrent with others (the so-called indicators) by using
```julia
LLR.indicators(A)
```
The result will be an item × item symmetric sparse matrix (see [SparseArrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) 
that contains all non-zero cooccurrence scores. Typically, you would pick the largest 50 or so items
in each row to be the indicators for an item.

## Example
For example, we can construct a matrix with some very common items and many items that are much less common
with 
```julia
using SparseArrays
p = collect(10 ./ (50 .+ (1:10000)))
ix = []; jx = []
for i in 1:10000
  j1 = findnz(sprand(Bool, 10000, p[i]))[1]
  n = length(j1)
  append!(ix, repeat([i], inner=n))
  append!(jx, j1)
end
A = sparse(ix, jx, 1)
```
The raw cooccurrence between items in this matrix is nasty because the first few rows have thousands of interactions
and that means that the cooccurrence requires many millions of computations. By bounding this, we get a very scalable
computation that is linear in execution time with the number of rows times `rowcut`. For this example, it takes about
20 seconds to find the indicators for `A` with `rowcut=100` on my laptop

