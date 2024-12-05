function [T, a, P, rho] = computeAtmosphericData(altitude)

    if altitude < 0
        altitude = 0;
    end
    
    %% Pressure and temperature 
    % Pressure calculation are
    % based on barometric formula
    % https://en.wikipedia.org/wiki/Barometric_formula
    
    H0_vec = [0 11 20 32 47 51 71]*1e3;                             % [m]
    P0_vec = [101325 22632.10 5474.89 868.02 110.91 66.94 3.96];    % [Pa]
    T0_vec = [288.15 216.65 216.65 228.65 270.65 270.65 214.65];    % [K]
    L_vec = [0.0065 0 -0.001 -0.0028 0 0.0028 0.002];               % [K/m]

    g0 = 9.80665;       % [m/s^2]
    M0 = 0.0289644;     % [kg/mol]
    R =  8.3144598;     % J/(mol*K)

    idx = sum(altitude >= H0_vec);
    if idx < 1
        idx = 1;
    end
    hb = H0_vec(idx);
    Pb = P0_vec(idx);
    Tb = T0_vec(idx);
    Lb = L_vec(idx);

    if Lb == 0
        P = Pb*exp((-g0*M0*(altitude-hb))/(R*Tb));
        T = Tb;
    else
        P = Pb * (1 - Lb/Tb * (altitude-hb))^((g0*M0)/(R*Lb));
        T = Tb - Lb*(altitude-hb);
    end
    if P < 0
        P = 0;
    end
    if T < 0 
        T = 0;
    end

    %% Density

    % Density data -> rho = rho_0 * exp(-(h-h0)/H)
    h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 ...
        180 200 250 300 350 400 450 500 600 700 800 900 1000]'*1e3;
    rho0_vect = [1.225 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3 ...
        3.206*1e-4 8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 ...
        2.438*1e-8 8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 ...
        2.789*1e-10 7.248*1e-11 2.418*1e-11 9.158*1e-12 3.725*1e-12 ...
        1.585*1e-12 6.967*1e-13 1.454*1e-13 3.614*1e-14 1.170*1e-14 ...
        5.245*1e-15 3.019*1e-15]';
    H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
        7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 ...
        53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]'*1e3;

    idx = sum(altitude >= h0_vect);
    h0 = h0_vect(idx);
    rho0 = rho0_vect(idx);
    H = H_vect(idx);
    
    rho = rho0 * exp(-(altitude-h0)/H);

    %% Air speed

    gamma = 1.4;
    R = P0_vec(1)/(rho0_vect(1)*T0_vec(1));
    a = sqrt(gamma*R*T);
end