module missile
using Printf

export Integrate
export propMass

propMass = 2000.0
mdotOn = propMass / 106.9
mdotOff = 0.0
inertMass = 1000.0
thrust = 1000.0

#tburn = 107.1
const vecSz = 7

mutable struct coeff 
  k1::Float64
  k2::Float64
  k3::Float64
  k4::Float64
end

Coeff() = coeff(0.0, 0.0, 0.0, 0.0)

#rkCoeffs = Array{coeff, 1}(undef, 6)
rkCoeffs = [Coeff() for _ in 1:7]

function ffunc(t::Float64, integ::Vector{Float64})
  gamma = atan(integ[6], integ[4])
  engOn = integ[7] > 0.0
  # str = @sprintf("%f %d", t, engOn ? 1 : 0)
  # println(str)
  if engOn
    g =  0.0
    mdot= -mdotOn
    acc = thrust / (inertMass + integ[7])
  else
    g =  -9.81
    mdot= mdotOff
    acc = 0.0
  end  
  deriv =[integ[4], integ[5], integ[6], acc*cos(gamma), 0.0, acc*sin(gamma) + g, mdot]
end

function Integrate(t0::Float64, step::Float64, integ::Vector{Float64})::Vector{Float64}
  hv = copy(integ)
  # step 1
  deriv=ffunc(t0, hv)
  for idx in 1:vecSz
      rkCoeffs[idx].k1 = deriv[idx]*step
      hv[idx] = integ[idx] + rkCoeffs[idx].k1 / 2.0
  end
  # step 2
  deriv = ffunc(t0 + step / 2.0, hv)
  for idx in 1:vecSz
    rkCoeffs[idx].k2 = deriv[idx]*step
    hv[idx] = integ[idx] + rkCoeffs[idx].k2 / 2.0
  end
  # step 3
  deriv = ffunc(t0 + step / 2.0, hv)
  for idx in 1:vecSz
    rkCoeffs[idx].k3 = deriv[idx]*step
    hv[idx] = integ[idx] + rkCoeffs[idx].k3
  end
  # step 4
  deriv = ffunc(t0 + step, hv)
  for idx in 1:vecSz
    rkCoeffs[idx].k4 = deriv[idx]*step
    integ[idx] = integ[idx] + (rkCoeffs[idx].k1 + 2.0 * rkCoeffs[idx].k2 + 2.0 * rkCoeffs[idx].k3 + rkCoeffs[idx].k4) / 6.0;
  end
  integ
end



end