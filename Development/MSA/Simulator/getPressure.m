function P = getPressure(h)
    % Based on barometric formula
    % https://en.wikipedia.org/wiki/Barometric_formula
    
    H0_vec = [0 11 20 32 47 51 71]*1e3;                             % [m]
    P0_vec = [101325 22632.10 5474.89 868.02 110.91 66.94 3.96];    % [Pa]
    T0_vec = [288.15 216.65 216.65 228.65 270.65 270.65 214.65];    % [K]
    L_vec = [0.0065 0 -0.001 -0.0028 0 0.0028 0.002];               % [K/m]

    g0 = 9.80665;       % [m/s^2]
    M0 = 0.0289644;     % [kg/mol]
    R =  8.3144598;     % J/(mol*K)

    idx = sum(h >= H0_vec);
    hb = H0_vec(idx);
    Pb = P0_vec(idx);
    Tb = T0_vec(idx);
    Lb = L_vec(idx);

    if Lb == 0
        P = Pb*exp((-g0*M0*(h-hb))/(R*Tb));
    else
        P = Pb * (1 - Lb/Tb * (h-hb))^((g0*M0)/(R*Lb));
    end
    if P < 0
        P = 0;
    end
end