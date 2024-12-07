function [X_COM_t,Jyaw_t,Jpitch_t,Jroll_t,GEOMETRY] = MASS_PROPERTIES_2(M_prop_t,GEOMETRY)


GEOMETRY.m_lox_2t = M_prop_t * GEOMETRY.OF2 / (1+GEOMETRY.OF2);%[kg] mass of lox
GEOMETRY.m_fuel_2t =  M_prop_t * 1  / (1+GEOMETRY.OF2);%[kg] mass of rp1

x_fuel_20 = (GEOMETRY.b1+GEOMETRY.b2+(GEOMETRY.b3/2));
x_fuel_2 = x_fuel_20 - (((abs(GEOMETRY.m_fuel_2-GEOMETRY.m_fuel_2t))/GEOMETRY.m_fuel_2)*x_fuel_20);

x_lox_20 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+(GEOMETRY.b4/2);
x_lox_2 = x_lox_20 - (((abs(GEOMETRY.m_lox_2-GEOMETRY.m_lox_2t))/GEOMETRY.m_lox_2)*x_lox_20);

x_b1 = GEOMETRY.b1*2/3; % fairing+pay+adapt
x_b2=GEOMETRY.b1+ GEOMETRY.b2/2; % forward skirt 1
x_b3 = GEOMETRY.b1+GEOMETRY.b2+(GEOMETRY.b3/2); % tank fuel 2
x_b_fuel_2 = x_fuel_2;
x_b4 =  GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+(GEOMETRY.b4/2); % tank lox 2
x_b_lox_2 = x_lox_2;
x_b5 = GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+(GEOMETRY.b5/2); %aft skirt 2

m_fuel_2t =  M_prop_t * GEOMETRY.OF2 / (1+GEOMETRY.OF2);%[kg] mass of lox;
m_lox_2t =  M_prop_t * 1  / (1+GEOMETRY.OF2);%[kg] mass of rp1;

if m_fuel_2t<=0

    x_b_fuel_2=0;
    m_fuel_2t=0;

end

if m_lox_2t<=0

    x_b_lox_2=0;
    m_lox_2t=0;

end


x_vec = [x_b1;x_b2;x_b3;x_b_fuel_2;x_b4;x_b_lox_2;x_b5];

m_vec = [GEOMETRY.m1;GEOMETRY.m2;GEOMETRY.m_tank_f_2;m_fuel_2t;GEOMETRY.m_tank_lox_2;m_lox_2t;GEOMETRY.m5];

X_COM_t = (dot(x_vec,m_vec)/(sum(m_vec)));
GEOMETRY.X_COM_t=X_COM_t;

xb1 = abs(x_b1-X_COM_t);
xb2 = abs(x_b2-X_COM_t);
xb3 = abs(x_b3-X_COM_t);
xb_fuel_2=abs(x_b_fuel_2-X_COM_t);
xb4 = abs(x_b4-X_COM_t);
xb_lox_2=abs(x_b_lox_2-X_COM_t);
xb5 = abs(x_b5-X_COM_t);


if m_fuel_2t>=GEOMETRY.m_fuel_1

    xb_fuel_2=0;
    m_fuel_2t=0;

end

if m_lox_2t>=GEOMETRY.m_lox_2

    xb_lox_2=0;
    m_lox_2t=0;

end


x_vec_J = [xb1;xb2;xb3;xb_fuel_2;xb4;xb_lox_2;xb5];


x_quad = x_vec_J.*x_vec_J;


J = dot(m_vec,x_quad);
Jyaw_t=J;
Jpitch_t=J;

R2 = GEOMETRY.Diam_2/2;

Rb1 = R2;
Rb2 = R2;
Rb3 = R2;
Rb_fuel_2=R2;
Rb4 = R2;
Rb_lox_2=R2;
Rb5 = R2;


if m_fuel_2t>=GEOMETRY.m_fuel_2

    x_b_fuel2=0;
    m_fuel_2t=0;
    Rb_fuel_2=0;

end

if m_lox_2t>=GEOMETRY.m_lox_2

    x_b_lox2=0;
    m_lox_2t=0;
Rb_lox_2=0;
end



r_vec = [Rb1;Rb2;Rb3;Rb_fuel_2;Rb4;Rb_lox_2;Rb5];

r_vec_quad = r_vec.*r_vec;

Jroll_t = dot(r_vec_quad,m_vec);

GEOMETRY.Jyaw_t= Jyaw_t;
GEOMETRY.Jpitch_t = Jpitch_t;
GEOMETRY.Jroll_t = Jroll_t;


end