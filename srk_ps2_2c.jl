#SRK - 2c

#Set-up the plots
using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("PyPlot")

using LinearAlgebra

#----------------------------------------------------------------------
#Set up Intial Parameters
#----------------------------------------------------------------------

#Set initial values for x
x0 = [0.0; #m1
      0.0; #m2
      0.0; #m3
      0.0; #p1
      0.0; #p2
      0.0; #p3
      ]

#Code the initial values for r
TX0_m1 = 0; #placeholder
TX0_m2 = 0; #placeholder
TX0_m3 = 0; #placeholder
r0 = [ TX0_m1; #m1
      TX0_m2; #m2
      TX_m1; #m3
      0.0; #p1
      0.0; #p2
      0.0; #p3
      ]

#-----------------------------------------------------------------
#Set-up of Parameters
#-----------------------------------------------------------------

#Input the known values from problem statement
Lx_char = 1000; #nt
Lt_char = 333; #AA
Lx_1 = 1200; #nt
Lx_2 = 2400; #nt
Lx_3 = 600; #nt
Ll_1 = Lx_1/3; #AA
Ll_2 = Lx_2/3; #AA
Ll_3 = Lx_3/3; #AA
doubletime = 30; #minutes
dil_rate = .693/doubletime; #min^-1
mass_cell_water = 0.7 #percent
copies_per_cell = 200; #200 copies

#Degradation rate of mRNA (kdeg_m)
#global half life of mRNA = 5 min https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111927&ver=2&trm=mrna+half+life+e+coli&org=
halflife = 5;  #min
#convert to seconds
halflife2 = .022*60; #min
#convert to degradation rate, assuming first order kinetics
kdeg_m = .693/halflife2 #min-1
#Degradation rate of protein (kdeg_p), bionumbers
halflife3 = 20*60 #min
kdeg_p = .693/halflife3 #min-1

#Concentration of gene (G) #umol/gDW
#vol_cell = 6.7*10.0^-10*10.0^-6; #L/cell (Bionumbers, see parameter estimates)
#https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108815&ver=1&trm=volume+of+e+coli+cell&org=
vol_cell = 9e-17 # L/cell, from PSET1 soln BIND: 114922
dens_cell = 1.7; #g/L (Bionumbers)
#https://bionumbers.hms.harvard.edu/bionumber.aspx?id=109836&ver=5&trm=glazyrina+cell+mass&org=
mol_gene = copies_per_cell / (6.02 * 10^23) #mol
dryweight_cell = vol_cell*dens_cell*(1-mass_cell_water)
G_gDW = mol_gene/dryweight_cell*10.0^6 #umol/gDW
G_mM = mol_gene*1000/(vol_cell*mass_cell_water) #mM

#kej for each mi
e_X = 42; # nt/sec, rate
ke_mchar = e_X/Lx_char*60 #min-1, characteristic length
ke_m1 = ke_mchar * Lx_char/Lx_1 #min^-1
ke_m2 = ke_mchar * Lx_char/Lx_2 #min^-1
ke_m3 = ke_mchar * Lx_char/Lx_3 #min^-1
#kej for each pi
e_L = 14.5*60 #AA/min
#https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100233&ver=8&trm=translation+rate+e+coli&org=
#Dalbow DG, Young R. Synthesis time of beta-galactosidase in Escherichia coli B/r as a function of growth rate. Biochem J. 1975 Jul150(1):13-20. abstract, p.15 table 2 and p.19 right column 3rd paragraph
ke_pchar = e_L/Lt_char #min^-1
ke_p1 = ke_pchar * Lt_char/Ll_1
ke_p2 = ke_pchar * Lt_char/Ll_2
ke_p3 = ke_pchar * Lt_char/Ll_3

#R_X and R_L
#RNAP concentration = 30 nM
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100194&ver=8&trm=rnap+e+coli+M&org=
R_X=30*10^-6; #mM
R_X_gDW = R_X*1000*(vol_cell*mass_cell_water)/dryweight_cell #umol/gDW
#Ribosomes per cell, https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108946&ver=4&trm=ribosomes+in+e+coli+concentration&org=
rib_per_cell = 27000 #ribosomes/cell
rib_per_cell_umol = rib_per_cell/(6.02*10^23)*10^6 #umol/cell
#Ribosome concentration
R_L = rib_per_cell/(6.02*10^23)/(vol_cell*mass_cell_water)*1000 #mM
R_L_gDW = rib_per_cell_umol/dryweight_cell #umol/gDW

#kI for m or p (Bionumbers)
#Closed to open complex = k2 = .024 s-1 (McClure Paper)
kI_m = 1/25*60; #min-1, for mRNA m
kI_p = 5 #min-1 for protein p (Bionumbers, https://bionumbers.hms.harvard.edu/bionumber.aspx?id=112001&ver=3&trm=initiation+rate+translation+e+coli&org=)
#Pai A, You L: Optimal tuning of bacterial sensing potential. Molecular systems biology 2009, 5(286) :286 doi: 10.1038/msb.2009.43. Supplementary Text 2 p.23 top paragraph

#Set up K_X and K_T
#K_X
#Find K_x,j = McClureSlope*kI (see paper), converted to units of [mM]
McClureSlope = 0.12*10^-3*1/60; #mM*min
K_X = McClureSlope*kI_m #mM
#K_T
McClureSlope_eq = 1.04*10^-3*1/60; #mM*min, #Assumption, that is the ribosome slope
K_T = kI_p*10 #*1000 (mutlipy by factor for troublshooting) #mM

#set up tau_m (for mRNA) and tau_p (for protein)
tau_m1 = ke_m1/kI_m
tau_m2 = ke_m2/kI_m
tau_m3 = ke_m3/kI_m
tau_p1 = ke_p1/kI_p
tau_p2 = ke_p2/kI_p
tau_p3 = ke_p3/kI_p

#-----------------------------------------------------------------
#Set-up promoter control rates
#-----------------------------------------------------------------

#Set-up subfunction for f fractions
n = 1.5 #assumption for all, from pset 1
w1_pset1 = 0.0; #for constant background translation, Assume none
w2_pset1 = 300; #for I on p1
w1_assume = 300; #for (p1 on p2) and (p1 on p3)
w2_assume = 100; #for p2 on p3
w3_assume = 100; #for p3 on p2
w4_assume = 100; #extra
kc_pset1 = 300; #mM, to use in f
#f as a function of some input
f(j) = (j^n)/(kc_pset1^n + j^n)

#-----------------------------------------------------------------
#Set-up For-Loop for Iteration through the Values
#-----------------------------------------------------------------

#Time ranges
tstart = 0; #min
tstep = 0.1; #min
tend = 360; #min
tSim = collect(tstart:tstep:tend)

#Empty array for outputs
x_total = zeros(length(tSim),6);

#Set-up initial

#-----------------------------------------------------------------
#Set-up Inducer Concentration
#-----------------------------------------------------------------

#Set up time dependence on I
if t<60 #before 60 min
 I = 10 #mM
else #after 60 min
 I = 0 #mM
end

#Make I constant (for troublshooting)
#I=10



#Set-up the different u------------------------------------
#-------------For normal (not broken) leave these uncommented
u_m1 = (w1_pset1 + w2_pset1*f(I))/(1+w1_pset1 + w2_pset1*f(I))
u_m2 = (w1_pset1 + w1_assume*f(p1) + w3_assume*f(p3))/(1+w1_pset1 + w1_assume*f(p1) + w3_assume*f(p3))
u_m3 = (w1_pset1 + w1_assume*f(p1) + w2_assume*f(p2))/(1+w1_pset1 + w1_assume*f(p1) + w2_assume*f(p2))
#--------------For broken leave these uncommented
#u_m1 = (w1_pset1 + w2_pset1*f(I))/(1+w1_pset1 + w2_pset1*f(I))
#u_m2 = (w1_pset1 + w1_assume*f(p1) + w3_assume*f(p3))/(1+w1_pset1 + w1_assume*f(p1) + w3_assume*f(p3))
#u_m3 = (w1_pset1 + w1_assume*f(p1))/(1+w1_pset1 + w1_assume*f(p1) + w2_assume*f(p2))

#for proteins assume translation happens at the kinetic limit
u_p1 = 1
u_p2 = 1
u_p3 = 1

#-----------------------------------------------------------------
#Set-up TX and TL rates
#-----------------------------------------------------------------

r_X1 = ke_m1*R_X*(G_mM/(tau_m1*K_X + (tau_m1+1)*G_mM))
r_X2 = ke_m2*R_X*(G_mM/(tau_m2*K_X + (tau_m2+1)*G_mM))
r_X3 = ke_m3*R_X*(G_mM/(tau_m3*K_X + (tau_m3+1)*G_mM))
r_L1 = ke_p1*R_L*(m1/(tau_p1*K_T + (tau_p1+1)*m1))
r_L2 = ke_p2*R_L*(m2/(tau_p2*K_T + (tau_p2+1)*m2))
r_L3 = ke_p3*R_L*(m3/(tau_p3*K_T + (tau_p3+1)*m3))

TX_1 = r_X1 * u_m1
TX_2 = r_X2 * u_m2
TX_3 = r_X3 * u_m3
TL_1 = r_L1 * u_p1
TL_2 = r_L2 * u_p2
TL_3 = r_L3 * u_p3



#----------------------------------------------------------------------
#Plotting
#----------------------------------------------------------------------

using PyPlot
figure(figsize=(4*10,3*10))
figure(1)
plot(t,p1_gDW,color="black")
plot(t,p2_gDW,color="blue")
plot(t,p3_gDW,color="red")
plot(t,I,color="green")
xlabel("time (min)")
ylabel("Concentration (umol/gDW)")
