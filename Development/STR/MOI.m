function [Jy,Jz]=MOI(h_it,M_it,M1,M2)


GEOMETRY = GEO(h_it,M_it,M1,M2);
L_tot = GEOMETRY.l_tot;
x_com = L_tot-GEOMETRY.x_com;
x_fuel_20 = (L_tot-GEOMETRY.b1-GEOMETRY.b2-(GEOMETRY.b3/2));
x_lox_20 = (L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-(GEOMETRY.b4/2));
x_fuel_10 = (L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-GEOMETRY.b4-GEOMETRY.b5-GEOMETRY.b6-(GEOMETRY.b7/2));
x_lox_10 = (L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-GEOMETRY.b4-GEOMETRY.b5-GEOMETRY.b6-GEOMETRY.b7-(GEOMETRY.b8/2)-GEOMETRY.b78);

x_b1 = abs((L_tot-GEOMETRY.b1*2/3) - x_com); % Fairing+Payload
x_b2 = abs((L_tot-GEOMETRY.b1-(GEOMETRY.b2/2)) - x_com); % forward skirt
x_b3 = abs((L_tot-GEOMETRY.b1-GEOMETRY.b2-(GEOMETRY.b3/2)) - x_com); % fuel tank2
x_b_fuel2 = abs(x_fuel_20 - x_com); % fuel 2
x_b4 =  abs((L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-(GEOMETRY.b4/2)) - x_com); % lox tank2
x_b_lox2 = abs(x_lox_20 - x_com); % lox 2
x_b5 = abs((L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-GEOMETRY.b4-(GEOMETRY.b5/2)) - x_com); %aft skirt 2
x_b6 = abs((L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-GEOMETRY.b4-GEOMETRY.b5-(GEOMETRY.b6/2)) - x_com); % interstage
x_b7 = abs((L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-GEOMETRY.b4-GEOMETRY.b5-GEOMETRY.b6-(GEOMETRY.b7/2)) - x_com); % fuel tank1
x_b_fuel1 = abs(x_fuel_10 - x_com); % fuel 1
x_b8 = abs((L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-GEOMETRY.b4-GEOMETRY.b5-GEOMETRY.b6-GEOMETRY.b7-(GEOMETRY.b8/2)-GEOMETRY.b78) - x_com); % lox tank1
x_b_lox1 =abs(x_lox_10 - x_com); % lox 1
x_b9 = abs((L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-GEOMETRY.b4-GEOMETRY.b5-GEOMETRY.b6-GEOMETRY.b7-GEOMETRY.b8-GEOMETRY.b78-(GEOMETRY.b9/2)) - x_com); % aft skirt 
x_b10 = abs((L_tot-GEOMETRY.b1-GEOMETRY.b2-GEOMETRY.b3-GEOMETRY.b34-GEOMETRY.b4-GEOMETRY.b5-GEOMETRY.b6-GEOMETRY.b7-GEOMETRY.b8-GEOMETRY.b78-GEOMETRY.b9-(GEOMETRY.b10/2)) - x_com); % engine 1

x_vec = [x_b1,x_b2,x_b3,x_b4,x_b5,x_b6,x_b7,x_b8,x_b9,x_b10];

x_quad = x_vec.*x_vec;

m_vec = [GEOMETRY.m1;GEOMETRY.m2;GEOMETRY.m3;GEOMETRY.m4;GEOMETRY.m5;GEOMETRY.m6;GEOMETRY.m7;GEOMETRY.m8;GEOMETRY.m9;GEOMETRY.m10];

J = dot(m_vec,x_quad);
Jy=J;
Jz=J;






end