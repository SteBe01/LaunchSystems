function [MAT] = material_selection(mat_id)


%this function gets the info about the selected material:
% mat_id = 1 for Ti, 
% 2 for Al 2XXX, 
% 3 for Steel, 
% 4 for Carbon Fiber,
% 5 for Al 7XXX
% 6 for Al-Li 2090

%switch
switch mat_id 
    case 1 % Ti6Al4V
        rho = 4500; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 110 * 1e9; %[Pa] young modulus
        sy = 900 * 1e6; %[Pa] tensile yield stress
        su = 950 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.34; %[-] Poisson's ratio
        Name = 'Ti6Al4V';
    case 2 % Al 2XXX
        rho = 2700; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 290 * 1e6; %[Pa] tensile yield stress
        su = 390 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
        Name = 'Al (2XXX)';
    case 3 % Steel
        rho = 7800; %[kg/m^3]
        t_min = 0.25 * 1e-3; %[m] minimum thickness for manufacturability
        E = 200 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = 420 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.27; %[-] Poisson's ratio
         Name = 'Steel';
    case 4 % Carbon fiber
        rho = 1800; %[kg/m^3]
        t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
        E = 250 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
        nu = 0.27; %[-] Poisson's ratio
        Name = 'Carbon Fiber';
    case 5 % Al 7XXX
        rho = 2750; %[kg/m^3]
        t_min = 1.06 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 500 * 1e6; %[Pa] tensile yield stress
        su = 510 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
         Name = 'Al (7XXX)';
    case 6 % AlLi (2090)
        rho = 2590; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability  
        E = 76 * 1e9; %[Pa] young modulus
        sy = 500 * 1e6; %[Pa] tensile yield stress
        su = 550 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.34; %[-] Poisson's ratio
        Name = 'Al-Li (2090)';
end

%recover material properties:
MAT.Name = Name;
MAT.ID = mat_id;%[-] ID of the material: 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX
MAT.rho = rho; %[kg/m^3] material density
MAT.t_min = t_min; %[m] material manufacturability minimum thickness
MAT.E = E; %[Pa] Young modulus
MAT.sy = sy; %[Pa] yelding stress
MAT.su = su; %[Pa] ultimate stress
MAT.nu = nu; %[-] Poisson's ratio



end