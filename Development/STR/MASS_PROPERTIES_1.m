function [X_COM_t,Jyaw_t,Jpitch_t,Jroll_t,GEOMETRY] = MASS_PROPERTIES_1(M_prop_t,GEOMETRY)


GEOMETRY.m_lox_1t = M_prop_t * GEOMETRY.OF1 / (1+GEOMETRY.OF1);%[kg] mass of lox
GEOMETRY.m_fuel_1t =  M_prop_t * 1  / (1+GEOMETRY.OF1);%[kg] mass of rp1

x_fuel_10 = (GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+GEOMETRY.b5+GEOMETRY.b6+(GEOMETRY.b7/2));
x_fuel_1 = x_fuel_10 - (((abs(GEOMETRY.m_fuel_1-GEOMETRY.m_fuel_1t))/GEOMETRY.m_fuel_1)*x_fuel_10);

x_lox_10 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+GEOMETRY.b5+GEOMETRY.b6+GEOMETRY.b7+(GEOMETRY.b8/2)+GEOMETRY.b78;
x_lox_1 = x_lox_10 - (((abs(GEOMETRY.m_lox_1-GEOMETRY.m_lox_1t))/GEOMETRY.m_lox_1)*x_lox_10);

R1 = GEOMETRY.Diam_1/2;
R2 = GEOMETRY.Diam_2/2;

if R1==R2

x_com_frustum = GEOMERTY.b6/2;

elseif R1~=R2

h=GEOMETRY.b6;
x_com_frustum = (h / 4) * ((R1^2 + 2*R1*R2 + 3*R2^2) / (R1^2 + R1*R2 + R2^2));

end

x_b1 = GEOMETRY.b1*2/3; % fairing+pay+adapt
x_b2=GEOMETRY.b1+ GEOMETRY.b2/2; % forward skirt 1
x_b3 = GEOMETRY.b1+GEOMETRY.b2+(GEOMETRY.b3/2); % tank fuel 2
%x_b_fuel_2 = GEOMETRY.b1+GEOMETRY.b2+((GEOMETRY.b3/2)); % fuel 2
x_b4 =  GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+(GEOMETRY.b4/2); % tank lox 2
%x_b_lox2 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+(GEOMETRY.b4/2); % lox 2
x_b5 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+(GEOMETRY.b5/2); %aft skirt 2
x_b6 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+GEOMETRY.b5+(x_com_frustum); % interstage
x_b7 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+GEOMETRY.b5+GEOMETRY.b6+(GEOMETRY.b7/2); % fuel tank1
x_b_fuel1 = x_fuel_1; % fuel 1
x_b8 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+GEOMETRY.b5+GEOMETRY.b6+GEOMETRY.b7+(GEOMETRY.b8/2)+GEOMETRY.b78; % lox tank1
x_b_lox1 =x_lox_1; % lox 1
x_b9 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+GEOMETRY.b5+GEOMETRY.b6+GEOMETRY.b7+GEOMETRY.b8+GEOMETRY.b78+(GEOMETRY.b9/2); % aft skirt 
x_b10 =GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+GEOMETRY.b5+GEOMETRY.b6+GEOMETRY.b7+GEOMETRY.b8+GEOMETRY.b78+GEOMETRY.b9+(GEOMETRY.b10/2); % engine 1

m_fuel_baseline =  GEOMETRY.m_prop_1 * 1  / (1+GEOMETRY.OF1);%[kg] mass of rp1
m_lox_baseline =  GEOMETRY.m_prop_1 * GEOMETRY.OF1 / (1+GEOMETRY.OF1);%[kg] mass of lox

m_fuel_1t =  M_prop_t * GEOMETRY.OF1 / (1+GEOMETRY.OF1);%[kg] mass of lox;
m_lox_1t =  M_prop_t * 1  / (1+GEOMETRY.OF1);%[kg] mass of rp1;

if M_prop_t>GEOMETRY.m_prop_1

    x_b_fuel1=x_fuel_10;
    m_fuel_1t=m_fuel_baseline;
       x_b_lox1=x_lox_10;
    m_lox_1t=m_lox_baseline;
 fprintf('Error, too much propellant1\n');
end


if m_fuel_1t<=0

    x_b_fuel1=0;
    m_fuel_1t=0;

end

if m_lox_1t<=0

    x_b_lox1=0;
    m_lox_1t=0;

end



x_vec = [x_b1;x_b2;x_b3;x_b4;x_b5;x_b6;x_b7;x_b_fuel1;x_b8;x_b_lox1;x_b9;x_b10];

m_vec = [GEOMETRY.m1;GEOMETRY.m2;GEOMETRY.m3;GEOMETRY.m4;GEOMETRY.m5;GEOMETRY.m6;GEOMETRY.m_tank_f_1;m_fuel_1t;GEOMETRY.m_tank_lox_1;m_lox_1t;GEOMETRY.m9;GEOMETRY.m10];


X_COM_t = (dot(x_vec,m_vec)/(sum(m_vec)));
GEOMETRY.X_COM_t=X_COM_t;

xb1 = abs(x_b1-X_COM_t);
xb2 = abs(x_b2-X_COM_t);
xb3 = abs(x_b3-X_COM_t);
xb4 = abs(x_b4-X_COM_t);
xb5 = abs(x_b5-X_COM_t);
xb6 = abs(x_b6-X_COM_t);
xb7 =abs(x_b7-X_COM_t);
xb_fuel1=abs(x_b_fuel1-X_COM_t);
xb8 =abs(x_b8-X_COM_t);
xb_lox1=abs(x_b_lox1-X_COM_t);
xb9 =abs(x_b9-X_COM_t);
xb10 = abs(x_b10-X_COM_t);

if m_fuel_1t>=GEOMETRY.m_fuel_1

    xb_fuel1=0;
    m_fuel_1t=0;

end

if m_lox_1t>=GEOMETRY.m_lox_1

    xb_lox1=0;
    m_lox_1t=0;

end


x_vec_J = [xb1;xb2;xb3;xb4;xb5;xb6;xb7;xb_fuel1;xb8;xb_lox1;xb9;xb10];


x_quad = x_vec_J.*x_vec_J;

J = dot(m_vec,x_quad);
Jyaw_t=J;
Jpitch_t=J;

Rb1 = R2;
Rb2 = R2;
Rb3 = R2;
Rb4 = R2;
Rb5 = R2;
Rb6 = (R1+R2)/2;
Rb7 = R1;
Rb_fuel_1=R1;
Rb8= R1;
Rb_lox_1=R1;
Rb9 = R1;
Rb10 = R1;

if m_fuel_1t>=GEOMETRY.m_fuel_1

    x_b_fuel1=0;
    m_fuel_1t=0;
    Rb_fuel_1=0;

end

if m_lox_1t>=GEOMETRY.m_lox_1

    x_b_lox1=0;
    m_lox_1t=0;
Rb_lox_1=0;
end



r_vec = [Rb1;Rb2;Rb3;Rb4;Rb5;Rb6;Rb7;Rb_fuel_1;Rb8;Rb_lox_1;Rb9;Rb10];

r_vec_quad = r_vec.*r_vec;

Jroll_t = dot(r_vec_quad,m_vec);

GEOMETRY.Jyaw_t= Jyaw_t;
GEOMETRY.Jpitch_t = Jpitch_t;
GEOMETRY.Jroll_t = Jroll_t;


end