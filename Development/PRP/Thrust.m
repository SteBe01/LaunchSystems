function [T, m_dot_ox, m_dot_fu] = Thrust(perc, z)

    % Calcola la spinta e il consumo di carburante in funzione di manetta e
    % altitudine interpolando i dati ottenuti da trottoling.m
    % 
    % INPUT:
    % perc  :   percentuale di manetta  [-]
    % z     :   altitudine              [m]
    % A_e   :   area di efflusso        [m^2]
    % 
    % OUTPUT:
    % T         :   spinta complessiva      [N]
    % m_dot_ox  :   portata di ossidante    [kg/s]
    % m_dot_fu  :   portata di fuel         [kg/s]

    data = importdata("data_throttle.mat");
    p0   = 101325;

    A_e = 0.0862;

    p_amb = p0*exp(-z/8600);
    p_exit = interp1(data.manetta, data.p_e, perc);
    T_static = A_e*(p_exit-p_amb);
    T = interp1(data.manetta, data.TT, perc) + T_static;
    m_dot_ox = interp1(data.manetta, data.m_dot_ox, perc);
    m_dot_fu = interp1(data.manetta, data.m_dot_fu, perc);

end