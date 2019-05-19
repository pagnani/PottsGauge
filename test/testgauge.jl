module TestGauge
using PottsGauge, Test


function gaugetests(q::Int,N::Int)

    J = rand(q,q,N,N)
    h = rand(q,N)
    J = 0.5(J + permutedims(J,[2,1,4,3] )) #symmetryzed

    Jtest,htest = gauge(J,h,ZeroSumGauge())
    @test PottsGauge.isgauge(Jtest,htest, ZeroSumGauge())
    μ,s=PottsGauge.testgauge(Jtest,htest, J,h)
    @test s/μ < 1e-8

    Jtest,htest = gauge(J,h,LatticeGas())
    @test PottsGauge.isgauge(Jtest,htest, LatticeGas())
    μ,s=PottsGauge.testgauge(Jtest,htest, J,h)
    @test s/μ < 1e-8

    xwt = rand(1:q,N)
    Jtest,htest = gauge(J,h,WildType(xwt))
    @test PottsGauge.isgauge(Jtest,htest, WildType(xwt))
    μ,s=PottsGauge.testgauge(Jtest,htest, J,h)
    @test s/μ < 1e-8
end

gaugetests(21,10)
gaugetests(4,3)
gaugetests(4,10)

printstyled("All TestGauge passed!\n",color=:green,bold=true)

end # end module TestGauge
