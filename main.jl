include("missile.jl")
import .missile
using Printf


function main(first)
  v0=1.0
  α = 50.0
  local integ =[0.0, 0.0, 0.0, v0*cosd(α), 0.0, v0*sind(α), missile.propMass]
  local t = 0.0
  local step = 0.01
  idx=0
  f = open("/temp/output.txt", "w")
  write(f,"t x y z vx vy vz pmass")
  write(f, "\n")
  prevX = 0.0
  prevZ = 0.0
  while true
    if first
      if integ[7]<missile.mdotOn*step
        currStep = integ[7]/missile.mdotOn
        @printf("step changed, t: %f, new step: %f [%f] \n", t, currStep, integ[7])
        first = false
      else
        currStep = step
      end
    else
      currStep = step
    end
    missile.Integrate(t, currStep, integ)
    t += currStep
    prtv = [[t] ; integ]
    fmtd = map(x ->@sprintf("%f", x), prtv)
    write(f, join(fmtd, " "))
    write(f, "\n")
    if (integ[3]<0.0 && integ[6]<0.0)
      p = -prevZ / (integ[3] * prevZ)
      q = 1.0 - p
      hit = q * prevX + p * integ[1]
      @printf("range: %f\n", hit)
      break
    end
    prevX=integ[1]
    prevZ=integ[3]
  end
  close(f)
  idx+=1
end

main(true)
main(false)
println("*** finished ****")
