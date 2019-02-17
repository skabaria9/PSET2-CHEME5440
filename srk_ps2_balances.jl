function srk_ps2_balances(t,x)

    #-----------------------------------------------------------------
    #Set-up of Species Vector
    #-----------------------------------------------------------------
    #Define x species vector
        m1 = x[1]
        m2 = x[2]
        m3 = x[3]
        p1 = x[4]
        p2 = x[5]
        p3 = x[6]

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
    halflife2 = 5; #min
    #convert to degradation rate, assuming first order kinetics
    kdeg_m = .693/halflife2 #min-1
    #Degradation rate of protein (kdeg_p)
    halflife3 = 20*60 #min
    kdeg_p = .693/halflife3 #min-1

    #Concentration of gene (G) #umol/gDW
    vol_cell = 6.7*10.0^-10 * 10.0^-6; #L/cell (Bionumbers, see parameter estimates)
    dens_cell = 1; #g/L (Assumption because mostly water)
    mol_gene = copies_per_cell / (6.02 * 10^23) #mol
    dryweight_cell = vol_cell*dens_cell*(1-mass_cell_water)
    G = mol_gene/dryweight_cell/10.0^6 #umol/gDW
    G_mM = mol_gene/(vol_cell*mass_cell_water)*1000 #mM

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
    #RNAP concentration
    # https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100194&ver=8&trm=rnap+e+coli+M&org=
    R_X=30*10^-6; #mM
    #Ribosomes per cell, https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108946&ver=4&trm=ribosomes+in+e+coli+concentration&org=
    rib_per_cell = 27000 #ribosomes/cell
    #Ribosome concentration
    R_L = rib_per_cell/(6.02*10^23)/(vol_cell*mass_cell_water)*1000 #mM

    #kI for m or p (Bionumbers)
    #Closed to open complex = k2 = .024 s-1 (McClure Paper)
    kI_m = 0.024*60; #min-1, for mRNA m
    kI_p = 5 #min-1 for protein p (Bionumbers, https://bionumbers.hms.harvard.edu/bionumber.aspx?id=112001&ver=3&trm=initiation+rate+translation+e+coli&org=)
    #Pai A, You L: Optimal tuning of bacterial sensing potential. Molecular systems biology 2009, 5(286) :286 doi: 10.1038/msb.2009.43. Supplementary Text 2 p.23 top paragraph

    #Set up K_X and K_T
    #K_X
    #Find K_x,j = McClureSlope*kI (see paper), converted to units of [mM]
    McClureSlope = 1.04*10^-3*1/60; #mM*min
    K_X = McClureSlope*kI_m #mM
    #K_T
    McClureSlope_eq = 1.04*10^-3*1/60; #mM*min, #Assumption
    K_T = McClureSlope_eq*kI_p #mM

    #set up tau_m (for mRNA) and tau_p (for protein)
    tau_m1 = ke_m1/kI_m
    tau_m2 = ke_m2/kI_m
    tau_m3 = ke_m3/kI_m
    tau_p1 = ke_p1/kI_p
    tau_p2 = ke_p2/kI_p
    tau_p3 = ke_p3/kI_p

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

    #-----------------------------------------------------------------
    #Set-up promoter control rates
    #-----------------------------------------------------------------

    #Set-up subfunction for f fractions
    n = 1.5 #assumption for all, from pset 1
    w1_pset1 = 0.0; #for constant background translation, Assume none
    w2_pset1 = 10; #for I
    w1_assume = 10;
    w2_assume = 10;
    w3_assume = 10;
    w4_assume = 10;
    kc_pset1 = 10;#mM
    #f as a function of some input
    f(j) = (j^n)/(kc_pset1^n + j^n)

    #Set-up the different u
    u_m1 = (w1_pset1 + w2_pset1*f(I))/(1+w1_pset1 + w2_pset1*f(I))
    u_m2 = (w1_pset1 + w1_assume*f(p1) + w3_assume*f(p3))/(1+w1_pset1 + w1_assume*f(p1) + w3_assume*f(p3))
    u_m3 = (w1_pset1 + w1_assume*f(p1) + w2_assume*f(p2))/(1+w1_pset1 + w1_assume*f(p1) + w2_assume*f(p2))
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

    #-----------------------------------------------------------------
    #Set-up Mass Balances(ODEs)
    #-----------------------------------------------------------------
    dxdt = similar(x)
    dxdt[1] = TX_1 - kdeg_m*m1 - dil_rate*m1 #m1
    dxdt[2] = TX_2 - kdeg_m*m2 - dil_rate*m2 #m2
    dxdt[3] = TX_3 - kdeg_m*m3 - dil_rate*m3#m3
    dxdt[4] = TL_1 - kdeg_p*p1 - dil_rate*p1#p1
    dxdt[5] = TL_2 - kdeg_p*p2 - dil_rate*p2#p2
    dxdt[6] = TL_3 - kdeg_p*p3 - dil_rate*p3 #p3
    dxdt

end
