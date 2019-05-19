# PottsGauge

Gauge transforms over the pair of Hamiltonian Potts parameters `(J,h)`. Assuming we have
a vector `x` of length `N`, and `x[i] ∈ 1:q` (`q` state Potts variable), a gauge transform
is a change `(J,h)→(J',h')` such that for all Potts configurations `x`
`E(x,J,h)-E(x,J',h')=const`.

The module defines  the gauge-transforms `gauge` of pair `J,h`. `J,h` are arrays of size `q×q×N×N, q×N` respectively. `J` should be symmetric: `J[a,b,i,j] == J[b,a,j,i]`. Allowed `gauge` are:  `ZeroSumGauge,LatticeGas,WildType`.

```
Usage:

    gauge(J,h,ZeroSumGauge())
    gauge(J,h,Wildtype(x0))
where in the last case `x0` is a `Vector{Int}` of size `N` whose values are in the interval `1,…,q`
```
