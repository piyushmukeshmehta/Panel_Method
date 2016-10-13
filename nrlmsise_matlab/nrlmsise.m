function [output] = nrlmsise(alt,lat,lon,year,doy,sec,lst,f107A,f107,ap_a)

input.alt = alt; % altitude in kilometers
input.g_lat = lat;
input.g_long = lon;
input.year = year;
input.doy = doy;
input.sec = sec;
input.lst = lst;
input.f107A = f107A;
input.f107 = f107;
% input.ap = ap;
input.ap_a = ap_a;

flags.switches(1)=1; % 0 - density in centimeters and grams
                     % 1 - density in meters and kilograms
for i=2:24
    flags.switches(i)=1;
end
flags.switches(10) = -1;
output = atmosnrlmsise00(input,flags); % centimeters and grams

% % matlab output
% t = 1.0e+003 * [1.0273    0.2212]; % kelvin
% rho = 1.0e+024 * [0.0000   0    6.6824    1.7927    0.0799    0.0000     0    0    0]; % kg/m^3
% 
% disp('temperature')
% output.t'
% t'
% disp('density')
% output.d'
% rho'

end