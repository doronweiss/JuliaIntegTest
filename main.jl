include("missile.jl")
import .missile
using Printf


function main()
  v0=100.0
  α = 50.0
  local integ =[0.0, 0.0, 0.0, v0*cosd(α), 0.0, v0*sind(α), missile.propMass]
  local t = 0.0
  local step = 0.01
  #typeof(integ)
  idx=0
  f = open("output.txt", "w") 
  while true
    missile.Integrate(t, step, false, integ)
    fmtd = map(x ->@sprintf("%f", x), integ)
    write(f, fmtd)
    t += step
    if (integ[3]<0.0 && integ[6]<0.0) || t>100
      break
    end
  end
  close(f)
  idx+=1
end

main()
