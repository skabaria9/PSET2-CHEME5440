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
    doubletime = 30*60; #seconds
    dilution = .693/doubletime; #s^-1
    mass_cell_water = 0.7 #percent
    copies_per_cell = 200; #200 copies

    #Calculate other constants needed for later analysis
    #Degradation rate of mRNA (kdeg_mRNA)
    #global half life of mRNA = 5 min https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111927&ver=2&trm=mrna+half+life+e+coli&org=
    halflife = 5;  #min
    #convert to seconds
    halflife2 = 5*60; #sec
    #convert to degradation rate, assuming first order kinetics
    kdeg_mRNA = .693/halflife2 #s-1

    #Concentration of gene (G) #umol/gDW
    vol_cell = 6.7*10.0^-10 * 10.0^-6; #L/cell
    dens_cell = 1; #g/L (Assumption because mostly water)
    mol_gene = copies_per_cell / (6.02 * 10^23) #mol
    dryweight_cell = vol_cell*dens_cell*(1-mass_cell_water)
    G = mol_gene/dryweight_cell/10.0^6 #umol/gDW

    #-----------------------------------------------------------------
    #Make I dependent on time
    #-----------------------------------------------------------------

    if t < 5
        I = 10 #mM, before 5 min
    else
        I = 0 #mM, after 5 min
    end
    #-----------------------------------------------------------------
    #Set-up of ODEs
    #-----------------------------------------------------------------



    #-----------------------------------------------------------------

end
