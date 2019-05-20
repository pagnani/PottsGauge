module PottsGauge

export gauge, ZeroSumGauge, LatticeGas, WildType

using Statistics

# start type definition

abstract type Gauge end

struct ZeroSumGauge <: Gauge end
struct LatticeGas <: Gauge end
struct WildType <: Gauge
    x0::Vector{Int}
end
struct ExternalGauge <: Gauge 
    U::Array{Float64,3}
    V::Array{Float64,3}
    C::Array{Float64,1}
end

"""
    gauge(J,h,gauge::T)

Return the gauge-transform `gauge` of pair `J,h`. `J,h` are arrays of size `q×q×N×N, q×N` respectively. `J` should be symmetric: `J[a,b,i,j] == J[b,a,j,i]`. Allowed `gauge` are:  `ZeroSumGauge,LatticeGas,WildType`.

Usage:

    gauge(J,h,ZeroSumGauge())
    gauge(J,h,Wildtype(x0))
where in the last case `x0` is a `Vector{Int}` of size `N` whose values are in the interval `1,…,q`
"""
function gauge(J::Array{F,4},h::Array{F,2},gauge::T) where F <: Real where T <: Gauge

    q,q,N,N = size(J)
    qh,Nh = size(h)
    Nh == N || error("Nh=$Nh != N=$N")
    qh == q || error("qh=$qh != q=$q")
    sum(abs2,(J-permutedims(J,[2,1,4,3]))) > 1e-10 && error("J should be symmetric: J - permutedims(J,[2,1,4,3] != 0)")

    JT = zero(J)

    U,V = zeros(F,q,N,N),zeros(F,q,N,N)
    UV!(U,V,J,gauge)
    for i in 1:N
        for j in 1:N
            if i != j
                for a in 1:q
                    for b in 1:q
                        JT[b,a,j,i] = J[b,a,j,i] + V[a,i,j] + U[b,i,j]
                    end
                end
            end
        end
    end

   hT = zero(h)
    for i in 1:N
        for a in 1:q
            hT[a,i] = h[a,i]
            for j = i+1:N
                hT[a,i] -= V[a,i,j]
            end
            for j = 1:i-1
                hT[a,i] -= U[a,j,i]
            end
        end
    end

    vecshift = if typeof(gauge) == WildType
        shift(J,h,gauge)
    else
        shift(hT,gauge)
    end
    for i in 1:N
        for a in 1:q
            hT[a,i] = hT[a,i] - vecshift[i]
        end
    end
    return JT,hT
end

"""
    function Energy(J,h,x)

Compute the energy with parameters `J::Array{T,4},h::Array{T,2}` of a  configuration `x::Array{Int,1}` of integer elements. `x[i] ∈ 1:q` where `q` is the number of letters (e.g. `q=21` in MSA with gaps);
"""
function Energy(J,h,x)

    q,q,N,N = size(J)
    qh,Nh = size(h)
    Nx = length(x)
    Nh == N || error("Nh=$Nh != N=$N")
    qh == q || error("qh=$qh != q=$q")
    Nx == Nh || error("Nx=$Nx != N=$Nh")
    all(z->1≤z≤q,x) || error("wrong alphabet")
    ene = 0.0
    @inbounds for i in 1:N-1
        for j in i+1:N
            ene -= J[x[j],x[i],j,i]
        end
        ene -= h[x[i],i]
    end
    ene -= h[x[N],N]
    return ene
end

"""
    function testgauge(J1,h1,J2,h2; nsample::Integer)
Test if parameters J1,h1 and J2,h2 are gauge related; Return mean and std of the energy difference induced by the two set of parameters J1,h1 and J2,h2 of `nsample` random configurations.
"""
function testgauge(J1,h1, J2,h2; nsample::Integer=100)
    q,q,N,N = size(J1)
    (size(J1) == size(J2)) || error("size J2 != size J1")
    (size(h1) == size(h2)) || error("size h2 != size h1")
    res = fill(0.0,nsample)
    for i=1:nsample
        conf = rand(1:q,N)
        E1 = Energy(J1,h1,conf)
        E2 = Energy(J2,h2,conf)
        res[i]= E2-E1
        # println("E1=$E1 E2=$E2\tE2-E1 =", E2-E1)
    end
    s = std(res)
    μ = mean(res)
    s > 1e-8 && warn("borked gauge? μ/std = ", μ/s)
    # warn("⟨E₂-E₁⟩ = ",μ,"\tstd = $s (nsample = $nsample)")
    return μ,s
end

function UV!(U::Array{T,3},V::Array{T,3},J::Array{T,4}, x::ZeroSumGauge) where T<:Real
    q,q,N,N = size(J)
    for i in 1:N
        for j in 1:N
            if i != j
                Jss = mean(J[:,:,i,j])
                V[:,i,j] = -vec(mean(J[:,:,i,j],dims=2)) .+ 0.5Jss
                U[:,i,j] = -vec(mean(J[:,:,i,j],dims=1)) .+ 0.5Jss
            end
        end
    end
end

function UV!(U::Array{T,3},V::Array{T,3},J::Array{T,4}, x::LatticeGas) where T<:Real
    q,q,N,N = size(J)
    for i in 1:N
        for j in 1:N
            if i != j
                for a in 1:q
                    V[a,i,j] = -J[a,q,i,j] + 0.5*J[q,q,i,j]
                    U[a,i,j] = -J[q,a,i,j] + 0.5*J[q,q,i,j]
                end
            end
        end
    end
end

function UV!(U::Array{T,3},V::Array{T,3},J::Array{T,4}, x::WildType) where T<:Real
    q,q,N,N = size(J)
    x0 = x.x0
    all(x -> 1 ≤ x ≤ q, x0) || error("gauge conf not compatible")
    for i in 1:N
        for j in 1:N
            if i != j
                for a in 1:q
                    V[a,i,j] = -J[a,x0[j],i,j] + 0.5 * J[x0[i],x0[j],i,j]
                    U[a,i,j] = -J[x0[i],a,i,j] + 0.5 * J[x0[i],x0[j],i,j]
                end
            end
        end
    end
end

function UV!(U::Array{T,3},V::Array{T,3},J::Array{T,4},x::ExternalGauge) where T<:Real
    V = x.V
    q,q,N,N = size(J)
    qU,NU,NU = size(U)
    qV,NV,NV = size(V)
    NU == N || error("size error on NU=$NU")
    NV == N || error("size error on NV=$NV")
    qU == q || error("size error on qU=$qU")
    qV == q || error("size error on qU=$qV")
    for i in 1:N
        for j in 1:N
            for a in 1:q
                U[a,i,j] = x.U[a,i,j]
                V[a,i,j] = x.V[a,i,j]
            end
        end
    end
end


shift(hT,::ZeroSumGauge) = ([ mean(hT[:,i]) for i in 1:size(hT,2)])
shift(hT,::LatticeGas) = ([ hT[end,i] for i in 1:size(hT,2)])
shift(ht,x::ExternalGauge) = (x.C)

function shift(J,h,x::WildType)
    N = length(x.x0)
    e0 = Energy(J,h,x.x0)
    [-e0/N for i in 1:size(h,2)]
end

function isgauge(J,h,::ZeroSumGauge)
    q,q,N,N = size(J)
    for i in 1:N
        for j in 1:N
            s = mean(J[:,:,i,j])
            if s > 1e-3
                @warn "s_{$i,$j} = $s ! not zsg"
                return false
            end
        end
    end
    return true
end

isgauge(J,h,g::WildType) = abs(Energy(J,h, g.x0)) > 1e-10 ? false : true

function isgauge(J,h,::LatticeGas)
    q,N=size(h)
    for i in 1:N-1
        for j in i+1:N
            Jb = view(J,:,:,i,j)
            srow = sum(abs2,Jb[:,end])
            scol = sum(abs2,Jb[end,:])
            scol + srow > 1e-10 && return false
        end
    end
    return true
end

end # end module GaugePotts
