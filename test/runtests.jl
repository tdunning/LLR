using Test
using LLR

# basic test values
@test LLR.g2test([1 0; 0 1]) ≈ 4 * log(2)
@test LLR.g2test([10 0; 0 10]) ≈ 40 * log(2)
@test LLR.g2test([10 20; 1 2]) ≈ 0 atol=1e-14
@test LLR.g2test([1 0; 10 0]) ≈ 0 
@test LLR.g2test([1 0; 0 20]) ≈ 8.040651442224128
@test LLR.g2test([1 0; 3 20]) ≈ 3.815168767010041

@test LLR.g2test([1 2 3 ; 4 5 6]) ≈ 0.29026921569868946

@test LLR.signedG2([1 0; 0 1]) ≈ sqrt(4 * log(2))
@test LLR.signedG2(1, 0, 0, 1) ≈ sqrt(4 * log(2))
@test LLR.signedG2(0, 1, 1, 0) ≈ -sqrt(4 * log(2))

@test LLR.signedG2.([1, 0, 10], [0, 1, 0], [0, 1, 0], [1, 0, 10]) ≈ [sqrt(4 * log(2)), -sqrt(4 * log(2)), sqrt(40 * log(2))]

# comparison of dictionaries
d = LLR.compare(Dict("a" => 1), Dict("b" => 1))
@test d["a"] ≈ sqrt(4 * log(2))
@test d["b"] ≈ -sqrt(4 * log(2))

d = LLR.compare(Dict("a" => 1, "b" => 1), Dict("b" => 10, "c" => 100))
@test d["a"] ≈ LLR.signedG2(1, 1, 0, 110)
@test d["b"] ≈ LLR.signedG2(1, 1, 10, 100)

# cooccurrence counting
A = [1 0 0 0 0; 0 2 0 1 0; 0 0 0 1 1; 0 1 0 1 1; 0 1 0 0 0]

# that matrix looks like this:
# 5×5 SparseMatrixCSC{Int64, Int64} with 9 stored entries:
#  1  .  .  .  .
#  .  2  .  1  .
#  .  .  .  1  1
#  .  1  .  1  1
#  .  1  .  .  .
# columns 2+4, 2+5 and 4+5 have cooccurrences. The 2x2 contingency tables
# for these colums are as follows:
# 2 1     1 2     2 1
# 1 1     1 1     0 2
# and the signed LLR scores are 0.37, -0.37 and 1.71 respectively

indicators = LLR.indicators(A)
@test indicators == indicators'
@test indicators[2,4] ≈ LLR.signedG2(2,1,1,1) 
@test indicators[2,5] ≈ LLR.signedG2(1,2,1,1) 
@test indicators[4,5] ≈ LLR.signedG2(2,1,0,2) 
