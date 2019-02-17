function Balances(t,x)

    #Define x species vector
    A = x[1]
    B = x[2]
    C = x[3]

    #Setup Parameters
    k1 = 1; #[] -> A
    k2 = 1; #A -> B
    k3 = 1; #B -> C
    kd_A = 0.5; #A -> []
    kd_B = 0.5; #B -> []
    kd_C = 0.5; #C -> []

    #Setup Mass Balances
    dxdt = similar(x)
    dxdt[1] = k1*(1-C/(0.1+C)) - k2*A - kd_A*A #A
    dxdt[2] = k2*A - k3*B - kd_B*B #B
    dxdt[3] = k3*B - kd_C*C #C
    #dxdt = [A;B;C]*[[-k2-kd_A 0.0 -0.1];[k2 -k3-kd_B 0.0];[0.0 k3 -kd_C]]
    dxdt
end
