function [V_fair] = V_fair_approx(D_fair_base_int,D_fair_base_ext,H_fair_base,D_fair_conetrap_int,D_fair_conetrap_ext,H_fair_conetrap,H_fair_nose)

V_fair.V_fair_base_int = (pi*(D_fair_base_int^2)/4)*H_fair_base; % [m^3]
V_fair.V_fair_base_ext = (pi*(D_fair_base_ext^2)/4)*H_fair_base; % [m^3]

V_fair.V_fair_conetrap_int = (1/3)*pi*H_fair_conetrap*( ((D_fair_base_int/2)^2) + ((D_fair_conetrap_int/2)^2) + ((D_fair_base_int/2)*(D_fair_conetrap_int/2)) ); % [m^3]
V_fair.V_fair_conetrap_ext = (1/3)*pi*H_fair_conetrap*( ((D_fair_base_ext/2)^2) + ((D_fair_conetrap_ext/2)^2) + ((D_fair_base_ext/2)*(D_fair_conetrap_ext/2)) ); % [m^3]

V_fair.V_fair_nose_int = (pi*(D_fair_conetrap_int^2)/4)*H_fair_nose; % [m^3]

V_fair.V_fair_tot = V_fair.V_fair_conetrap_ext + V_fair.V_fair_base_ext + V_fair.V_fair_nose_int;
V_fair.V_fair_tot_int = V_fair.V_fair_conetrap_int + V_fair.V_fair_base_int + V_fair.V_fair_nose_int;

V_fair.thickness = (D_fair_base_ext - D_fair_base_int)/2 ;



end