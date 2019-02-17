include("balances_ex.jl")

using Pkg
Pkg.add("ODE")
Pkg.add("PyPlot")

using ODE

#Setup Time Vector
tStart = 0.0
tStep = 0.1
tStop = 10.0
tSim = collect(tStart:tStep:tStop)

#Setup initial conditions
x0 = [0.0; #A
      0.0; #B
      0.0; #C
      ]

f(t,x) = Balances(t,x)
t,X = ode23s(f,x0,tSim; points=:specified)

A = [a[1] for a in X]
B = [a[2] for a in X]
C = [a[3] for a in X]

using PyPlot
#figure(figsize=(4,3))
figure(1)
plot(t,A,color="black")
plot(t,B,color="blue")
plot(t,C,color="red")
xlabel("time (h)")
ylabel("Concentration (mM)")
axis([0, 10, 0, 1])
tight_layout()
