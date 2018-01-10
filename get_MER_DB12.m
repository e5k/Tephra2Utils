% TError package
% This file was written by W. Degruyter and C. Bonadonna:
% Degruyter, W., & Bonadonna, C. (2012). Improving on mass flow rate estimates of volcanic eruptions. Geophys Res Lett, 39(16). doi:10.1029/2012GL052566
% Get MER from Degruyter and Bonadonna
% H     =   Height above the vent (m)
% Vmax  =   Wind at the tropopause
% H1    =   Tropopause height above the vent (m)
function Mdot = get_MER_DB12(H, Vmax, H1)


%constants 
g   = 9.81;         % gravitational acceleration (m s^-2) 22
z_1 = 2.8;          % maximum non-dimensional height (Morton et al. 23 1956) 24
R_d = 287;          % specific gas constant of dry air (J kg^-1 K^-25 1)1232 26
C_d = 998;          % specific heat capacity at constant pressure of dry 27 air (J kg^-1 K^-1) 28
C_s = 1250;         % specific heat capacity at constant pressure of 29 solids (J kg^-1 K^-1) 30
theta_0 = 1300;     % initial plume temperature (K) 
% entrainment 
alpha   = 0.1;      % radial entrainment coefficient 35
beta    = 0.5;      % wind entrainment coefficient 36

% height of the plume above vent (m)
dummyH = H; %0:10:40000;

% atmosphere temperature profile (Woods, 1988) 41
theta_a0 = 288;     % atmopshere temperature at the vent (K) 42
P_0      = 101325;  % atmopshere pressure at the vent (Pa) 43
rho_a0   = P_0/(R_d*theta_a0); % reference density atmosphere (kg m^-3) 44
%H1 = 12000;         % height of the tropopause above the vent (m) 45
H2 = 20000;         % height of the stratosphere above the vent (m) 46
tempGrad_1 = -6.5/1000; % temperature gradient in the troposphere (K m^-47 1) 48
tempGrad_2 = 0;     % temperature gradient between troposphere and startosphere (K m^-1) 50
tempGrad_3 = 2/1000;% temperature gradient in the stratosphere (K m^-51 1) 52

% reduced gravity (m s^-2) 54
gprime = g*(C_s*theta_0-C_d*theta_a0)/(C_d*theta_a0);
% average square buoyancy frequency Nbar^2 = Gbar across height of 
% the plume (s^-2) 
G1 = g^2/(C_d*theta_a0)*(1+C_d/g*tempGrad_1); 
G2 = g^2/(C_d*theta_a0)*(1+C_d/g*tempGrad_2); 
G3 = g^2/(C_d*theta_a0)*(1+C_d/g*tempGrad_3); 
Gbar = G1.*ones(size(dummyH)); 
Gbar(dummyH>H1) = (G1.*H1 + G2.*(dummyH(dummyH>H1)-H1))./dummyH(dummyH>H1); 
Gbar(dummyH>H2) = (G1.*H1 + G2.*(H2-H1) + G3.*(dummyH(dummyH>H2)-H2))./dummyH(dummyH>H2); 
Nbar = Gbar.^(1/2); 
% atmosphere wind profile (Bonadonna and Phillips, 2003) 71
%Vmax = 0; % maximum wind speed at the tropopause (m/s), change to e.g. 30 to see effect of wind
% average wind speed across height of the plume (m/s) 
Vbar = Vmax.*dummyH./H1./2; 
Vbar(dummyH>H1) = 1./dummyH(dummyH>H1).*(Vmax.*H1./2 + Vmax.*(dummyH(dummyH>H1)-H1) -0.9.*Vmax./(H2-H1).*(dummyH(dummyH>H1)-H1).^2./2); 
Vbar(dummyH>H2) = 1./dummyH(dummyH>H2).*(Vmax.*H1./2 + 0.55.*Vmax.*(H2-H1) + 0.1.*Vmax.*(dummyH(dummyH>H2)-H2)); 
% equation (6) in manuscript 
Mdot = pi*rho_a0/gprime*((2^(5/2)*alpha^2.*Nbar.^3./z_1.^4).*dummyH.^4 +(beta^2.*Nbar.^2.*Vbar/6).*dummyH.^3); 
