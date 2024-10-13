module missile

export propMass, tburn
export Integrate

propMass = 1000.0
inertMass = 1000.0
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

function ffunc(t::Float64, integ::Vector{Float64})::Vector{Float64}
  gamma = atan(integ[6], integ[4])
  engOn = integ[7] > 0.0
  if engOn
    g =  0.0
    thrust = 100.0
    mdot= -1000.0 / 106.9
    acc = thrust / (inertMass + integ[7])
  else
    g =  -9.81
    thrust = 0.0
    mdot= 0.0
    acc = 0.0
  end  
  deriv =[integ[4], integ[5], integ[6], acc*cos(gamma), 0.0, acc*sin(gamma) + g, mdot]
end

function Integrate(t0::Float64, step::Float64, exact::Bool, integ::Vector{Float64})
  hv = integ
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
    integ = integ[idx] + (rkCoeffs[idx].k1 + 2.0 * rkCoeffs[idx].k2 + 2.0 * rkCoeffs[idx].k3 + rkCoeffs[idx].k4) / 6.0;
  end
end



end