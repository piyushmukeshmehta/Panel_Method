% -------------------------------------------------------------------------
% Test nrlmsise-00 atmospheric model
% -------------------------------------------------------------------------
function test_nrlmsise00()

% comparison with the output of nrlmsise-00_test.c
% test1()

% comparison with the output of the MATLAB function atmosnrlmsise00
test2() 

% comparison with exercise that use the MATLAB function atmosnrlmsise00
% test3() % non torna il density plot

% Density/altitude profile
% test4()

% comparison with the Harris-Priester Model
%test5()

end
% -------------------------------------------------------------------------
function test1()
% comparison with the output of nrlmsise-00_test.c

% input values
for i=1:7
    aph.a(i)=100;
end
flags.switches(1)=0;
for i=2:24
    flags.switches(i)=1;
end
for i=1:17
    input(i).doy=172;
    input(i).year=0; % without effect
    input(i).sec=29000;
    input(i).alt=400;
    input(i).g_lat=60;
    input(i).g_long=-70;
    input(i).lst=16;
    input(i).f107A=150;
    input(i).f107=150;
    input(i).ap=4;
end
input(2).doy=81;
input(3).sec=75000;
input(3).alt=1000;
input(4).alt=100;
input(11).alt=0;
input(12).alt=10;
input(13).alt=30;
input(14).alt=50;
input(15).alt=70;
input(17).alt=100;
input(5).g_lat=0;
input(6).g_long=0;
input(7).lst=4;
input(8).f107A=70;
input(9).f107=180;
input(10).ap=40;
input(16).ap_a=aph;
input(17).ap_a=aph;
%evaluate 0 to 14
for i=1:15
    %gtd7(&input(i), &flags, &output(i));
    output(i) = atmosnrlmsise00(input(i),flags);
    %outpu(i) = tmp;
end
% evaluate 15 and 16
flags.switches(9)=-1;
for i=16:17
    %gtd7(&input(i), &flags, &output(i));
    output(i) = atmosnrlmsise00(input(i),flags);
end

% output type 1
for i=1:17
    fprintf('\n');
    for j=1:9
        fprintf('%E ',output(i).d(j));
    end
    fprintf('%E ',output(i).t(1));
    fprintf('%E \n',output(i).t(2));
    % DL omitted
end

% output type 2
for i=1:3
    fprintf('\n');
    fprintf('\nDAY   ');
    for j=1:5
        fprintf('         %3i',input((i-1)*5+j).doy);
    end
    fprintf('\nUT    ');
    for j=1:5
        fprintf('       %5.0f',input((i-1)*5+j).sec);
    end
    fprintf('\nALT   ');
    for j=1:5
        fprintf('        %4.0f',input((i-1)*5+j).alt);
    end
    fprintf('\nLAT   ');
    for j=1:5
        fprintf('         %3.0f',input((i-1)*5+j).g_lat);
    end
    fprintf('\nLONG  ');
    for j=1:5
        fprintf('         %3.0f',input((i-1)*5+j).g_long);
    end
    fprintf('\nLST   ');
    for j=1:5
        fprintf('       %5.0f',input((i-1)*5+j).lst);
    end
    fprintf('\nF107A ');
    for j=1:5
        fprintf('         %3.0f',input((i-1)*5+j).f107A);
    end
    fprintf('\nF107  ');
    for j=1:5
        fprintf('         %3.0f',input((i-1)*5+j).f107);
    end
    fprintf('\n\n');
    fprintf('\nTINF  ');
    for j=1:5
        fprintf('     %7.2f',output((i-1)*5+j).t(1));
    end
    fprintf('\nTG    ');
    for j=1:5
        fprintf('     %7.2f',output((i-1)*5+j).t(2));
    end
    fprintf('\nHE    ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(1));
    end
    fprintf('\nO     ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(2));
    end
    fprintf('\nN2    ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(3));
    end
    fprintf('\nO2    ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(4));
    end
    fprintf('\nAR    ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(5));
    end
    fprintf('\nH     ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(7));
    end
    fprintf('\nN     ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(8));
    end
    fprintf('\nANM 0 ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(9));
    end
    fprintf('\nRHO   ');
    for j=1:5
        fprintf('   %1.3e',output((i-1)*5+j).d(6));
    end
    fprintf('\n');
end
fprintf('\n');

end
% -------------------------------------------------------------------------
function test2()
% TEST #2
% from http://uk.mathworks.com/help/aerotbx/ug/atmosnrlmsise00.html

input.alt = 10; % altitude in kilometers
input.g_lat=45;
input.g_long=-50;
input.year=2007;
input.doy=4;
input.sec=0;
flags.switches(1)=1; % 0 - density in centimeters and grams
                     % 1 - density in meters and kilograms
for i=2:24
    flags.switches(i)=1;
end
output = atmosnrlmsise00(input,flags); % centimeters and grams

% matlab output
t = 1.0e+003 * [1.0273    0.2212]; % kelvin
rho = 1.0e+024 * [0.0000   0    6.6824    1.7927    0.0799    0.0000     0    0    0]; % kg/m^3

disp('temperature')
output.t'
t'
disp('density')
output.d'
rho'

end
% -------------------------------------------------------------------------
function test3()
% TEST #3
% from Spacecraft Interactions Project, 28 August 2013

flags.switches(1)=0;  % centimeters and grams
for i=2:24
    flags.switches(i)=1;
end

% Ex1
% Using the MSISE-00 model, plot the total mass density of the upper atmosphere from 200 to
% 1000 km altitude for F10.7 = 140 sfu and Ap = 15 nT. Plot the same parameter for F10.7 =
% 220 sfu and Ap = 15 nT. How does the total mass density as a function of altitude change
% with the F10.7 flux?

input.g_lat=0;
input.g_long=0;
input.year=2012;
input.doy=1;
input.sec=0;
input.ap=15;

alt_max = 1000;
alt_min = 200;
step = 1;
alt = [alt_min:step:alt_max]; % altitude in kilometers

f107 = [140 220];

figure(1)
for j=1:2
    input.f107=f107(j);
    for i=1:length(alt)
        input.alt = alt(i);
        output = atmosnrlmsise00(input,flags);
        dens(i) = output.d(6);
    end
    %semilogy(alt,dens*1e-3/1e-6,'-')   % kg m-3
    semilogy(alt,dens,'-')              % g cm-3
    hold on
end
xlabel('Altitude [km]')
ylabel('Total mass density [g/cm^3]')
title('Mass Density Upper Atmosphere')
legend('140 sfu','220 sfu')

clear input

% Ex2
% I also wanted to see if NRLMSISE-00 would portray the diurnal bulge of the atmosphere.
% Figure 3 below is a plot of latitude and longitude, 350 km at 0000 UTC. Since it’s
% midnight UTC, I’d expect to see a decrease in atmospheric density at the prime meridian,
% which Figure 3 clearly shows.

%Contour plot
%Atmos Density
%h=350km
%F10.7/a = 140 sfu, Ap =15
%2013 Jan 1

format long g
input.alt = 350;
input.year = 2013;
input.doy = 182;
input.sec = 64800;
input.f107 = 140;
input.ap = 15;

i=1;
for long=-180:180
j=1;
input.g_long = long;
for lat = -90:90
    input.g_lat = lat;
output = atmosnrlmsise00(input,flags);
T = output.t;
den(i,j)=output.d(6)*1e3;   % kg m-3
j=j+1;
end
i=i+1;
end
lat=-90:90;
long=-180:180;

figure(2)

[LA,LO]=meshgrid(lat,long);
pcolor(LO,LA,den)
shading interp
xlabel('Longitude, deg')
ylabel('Latidude, deg')
title('Atmosphere Density, 1 Jul 2013 18:00:00, 350 km, F10.7=140, Ap=15')

% worldmap world
% load coast
% [latcells, loncells] = polysplit(lat, long);
% numel(latcells)
% plotm(lat, long)
end
% -------------------------------------------------------------------------
function test4()
% Density/Altitude profile

input.g_lat=45;
input.g_long=-50;
input.year=2007;
input.doy=4;
input.sec=0;
flags.switches(1)=0;  % centimeters and grams
for i=2:24
    flags.switches(i)=1;
end

alt_max = 1000;
alt_min = 0;
step = 10;
alt = [alt_min:step:alt_max]; % altitude in kilometers

for i=1:length(alt)
    input.alt = alt(i);
    output = atmosnrlmsise00(input,flags);
    dens(i) = output.d(6);
end

figure(100)
plot(alt,dens,'.-')

end


function test5()
% Comparison with the Harris-Priester model

addpath('../')
addpath(genpath('../astrotoolbox'))

muE = astroConstants(13);
tsf = 86400;

R_equ   =   6378.137;         % Radius Earth [km]; WGS-84
r_Sun = [149600000 0 0]; % km
% r0 = R_equ+150;
% y = [1 0 0 ] * r0;
% dy = [0 1 0] * sqrt(muE/r0);
U1 = [1 0 0; 0 1 0; 0 0 1];
%t = mjd2mjd2000(55600)
t = 6.688382900299999e+03

Am = 0.1;      % [m^2/kg]
Cr =-1.21;      % abs(Cr) > 1, Cr <0
Cd = 2.2;       % 1.5 < Cd < 3.0 

% Unit conversion
Am = Am/1e6;  % [km^2/kg]

load('Ephemeris_Gaiar1.mat')
ephemeris = readEphHorizons(t,ephTable(1).ephemeris);
r_pE = ephemeris(1:3);       %planet distance from the Earth
r_pE = transpose(r_pE);
r_Sun = r_pE;

count = 1;
for i=150:10:1000
    
    alt(count) = i;
    
    r0 = R_equ + i;
    y = [1 0 0 ] * r0;
    dy = [0 1 0] * sqrt(muE/r0);

    % Drag from Harris-Priester model
    da(count,:) = drag_force(r_Sun,y(1:3),dy(1:3),U1',Am,1,Cd);
    norm_da(count) = norm(da(count,:));
    
    % Drag from NRLMSISE-00
    p.am = Am; p.Cd = Cd; p.f107 = 150;
    [da2(count,:) dens_gcm3(count)] = drag_nrlmsise00(y(1:3),dy(1:3),U1',p,t-2000);
    norm_da2(count) = norm(da2(count,:));

    date = mjd20002date(t-2000); % date = [Y, M, D, hrs, mn, sec];
    doy    = ymd2doy(date(1),date(2),date(3));
    sec    = hms2sec(date(4),date(5),date(6));
    date1 = [doy sec];
    
    f107 = p.f107;
    da3(count,:) = drag_force_nrlmsise00( y(1:3), dy(1:3), U1', Am, 1, Cd, date1, f107);
    norm_da3(count) = norm(da3(count,:));
    count = count+1;
end

figure(1)
plot(alt,norm_da,'.-b',alt,norm_da2,'.-r',alt,norm_da3,'.-k')
ylabel(strcat('|a|'))
xlabel('altitude [km]')
legend('HP','NRLMSISE-00','NRLMSISE-00/drag mex')

figure(2)
for i=1:3
    subplot(1,3,i)
    plot(alt,da(:,i),'.b',alt,da2(:,i),'.r',alt,da3(:,i),'.k')
    ylabel(strcat('a[',num2str(i),']'))
    xlabel('altitude [km]')
    legend('HP','NRLMSISE-00','NRLMSISE-00/drag mex')
end

figure(3)
plot(alt,dens_gcm3,'.-r')

end
