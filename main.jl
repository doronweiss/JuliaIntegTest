include("missile.jl")
import .missile

v0=100.0
α = 50.0
integ =[0.0, 0.0, 0.0, v0*cosd(α), 0.0, v0*sind(α), missile.propMass]
#typeof(integ)
while true
  v = missile.Integrate(0.0, 0.01, false, integ)
  integ = v
  if integ[3]<0.0 && integ[6]<0.0
    break
  end
end