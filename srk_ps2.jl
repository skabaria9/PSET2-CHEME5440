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

m1 = [a[1] for a in X] #mM
m2 = [a[2] for a in X] #mM
m3 = [a[3] for a in X] #mM
p1a = [a[4] for a in X] #mM
p2a = [a[5] for a in X] #mM
p3a = [a[6] for a in X] #mM

#Convert protein and mRNA concentration from mM to gDW


using PyPlot
figure(figsize=(4,3))
figure(1)
plot(t,p1a,color="black")
plot(t,p2a,color="blue")
plot(t,p3a,color="red")
xlabel("time (min)")
ylabel("Concentration (mM)")
#axis([0, 10, 0, 1])
#tight_layout()
