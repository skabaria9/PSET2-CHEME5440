include("srk_ps2_balances.jl")

using Pkg
Pkg.add("ODE")
Pkg.add("PyPlot")

using ODE

#Setup Time Vector
tStart = 0.0 #min
tStep = 1.0 #min
tStop = 500.0 #min
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
I = zeros(length(t)) #set default all zeros in I
for i in 1:length(t)
      if t[i] < 60
            print(t[i])
            I[i] = 1
      else
            I[i] = 0
      end
end

#Convert protein and mRNA concentration from mM to gDW
#set-up conversion factor formula
#vol_cell = 6.7*10.0^-10*10.0^-6; #L/cell (Bionumbers, see parameter estimates)
#https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108815&ver=1&trm=volume+of+e+coli+cell&org=
vol_cell = 9e-17 # L/cell, from PSET1 soln BIND: 114922
fraction_cell_water = 0.7 #percent
dens_cell = 1.7
dryweight_cell = vol_cell*dens_cell*(1-mass_cell_water)
mass_of_single_cell = 2.8e-13                   # g/cell BIND:103904
dryweight_cell2 = (1-fraction_cell_water)*mass_of_single_cell;
convertfactor(mM) = (mM*1000)*(vol_cell)/dryweight_cell2 #umol/gDW
#convert all values
m1_gDW = convertfactor(m1)
m2_gDW = convertfactor(m2)
m3_gDW = convertfactor(m3)
p1_gDW = convertfactor(p1a)
p2_gDW = convertfactor(p2a)
p3_gDW = convertfactor(p3a)

using PyPlot
figure(figsize=(4*10,3*10))
figure(1)
plot(t,p1_gDW,color="black")
plot(t,p2_gDW,color="blue")
plot(t,p3_gDW,color="red")
plot(t,I,color="green")
xlabel("time (min)")
ylabel("Concentration (umol/gDW)")


figure(figsize=(4*10,3*10))
figure(2)
plot(t,p1a,color="black")
plot(t,p2a,color="blue")
plot(t,p3a,color="red")
xlabel("time (min)")
ylabel("Concentration (mM)")
#axis([0, 10, 0, 1])
#tight_layout(
