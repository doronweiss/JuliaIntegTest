include("missile.jl")
import .missile
using Printf


function main()
  v0=100.0
  α = 50.0
  local integ =[0.0, 0.0, 0.0, v0*cosd(α), 0.0, v0*sind(α), missile.propMass]
  local t = 0.0
  local step = 0.001
  idx=0
  f = open("output.txt", "w")
  write(f,"t x y z vx vy vz pmass")
  write(f, "\n")
  while true
    missile.Integrate(t, step, integ)
    t += step
    prtv = [[t] ; integ]
    fmtd = map(x ->@sprintf("%f", x), prtv)
    write(f, join(fmtd, " "))
    write(f, "\n")
    if (integ[3]<0.0 && integ[6]<0.0) || t>100
      break
    end
  end
  close(f)
  idx+=1
end

main()
println("*** finished ****")
