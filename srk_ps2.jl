include("srk_ps2_balances.jl")

using Pkg
Pkg.add("ODE")
Pkg.add("PyPlot")

using ODE

#Setup Time Vector
tStart = 0.0 #min
tStep = 1.0 #min
tStop = 360.0 #min
tSim = collect(tStart:tStep:tStop)

#Setup initial conditions
x0 = [0.0; #m1
      0.0; #m2
      0.0; #m3
      0.0; #p1
      0.0; #p2
      0.0; #p3
      ]

f(t,x) = srk_ps2_balances(t,x)
t,X = ode23s(f,x0,tSim; points=:specified)

#Print lines for troubleshooting the code
#print(t)
#print(X)

m1 = [a[1] for a in X]
m2 = [a[2] for a in X]
m3 = [a[3] for a in X]
p1 = [a[4] for a in X]
p2 = [a[5] for a in X]
p3 = [a[6] for a in X]

using PyPlot
figure(figsize=(4,3))
figure(1)
plot(t,p1,color="black")
plot(t,p2,color="blue")
plot(t,p3,color="red")
xlabel("time (min)")
ylabel("Concentration (mM)")
#axis([0, 10, 0, 1])
#tight_layout()
