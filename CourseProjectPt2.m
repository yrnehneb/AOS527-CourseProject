% Benjamin Henry
% AOS527: Atmospheric Radiative Transfer
% Course Project (Part 2: Pre-Project)
%
% This script tackles both of the pre-project tasks.

%Pre-Project Task #1: Calculating Fluxes
bands = [0 160 560 630 700 800 900 990 1070 1200 1400 ...
    2200 2250 2400 3000];  %Bands in Wavenumber limits (cm^-1)
bands = bands*100; %Conversion from cm^-1 to m^-1
T_CIRC2 = 307.74; %Temperature in Kelvin
T_CIRC4 = 266.55; %Temperature in Kelvin
Emiss = 1.0;      %Surface Emissivity

%Create the Planck Function 
h = 6.63*10.^(-34); %Plank Constant (J*s)
c = 3*10.^8; %Speed of Light (m/s)
k = 1.38*10.^(-23); %Boltzmann Constant (m^2*kg/s^2*K)
sig = 5.67*10.^(-8); %Stefan Boltzman Constant
syms B_CIRC2(v);syms B_CIRC4(v) 
B_CIRC2(v) = 2*h*(c.^2)*(v.^3)*(1/(exp(h*c*v/(k*T_CIRC2))-1));
B_CIRC4(v) = 2*h*(c.^2)*(v.^3)*(1/(exp(h*c*v/(k*T_CIRC4))-1));
B_CIRC2 = matlabFunction(B_CIRC2);
B_CIRC4 = matlabFunction(B_CIRC4);


%Integrating Planck Function to Find Intensity
%CIRC2 Profile
for i = 1:length(bands)-1
    band_flux_CIRC2(i) = integral(B_CIRC2, bands(i), bands(i+1));
end
tot_band_flux_CIRC2 = sum(band_flux_CIRC2)*pi; %Converting from Intensity to Flux
tot_flux_CIRC2 = Emiss*sig*T_CIRC2.^4;

%CIRC4 Profile
for i = 1:length(bands)-1
    band_flux_CIRC4(i) = integral(B_CIRC4, bands(i), bands(i+1));
end
tot_band_flux_CIRC4 = sum(band_flux_CIRC4)*pi; %Converting from Intensity to Flux
tot_flux_CIRC4 = Emiss*sig*T_CIRC4.^4;

%Print Results
%CIRC2 Results
sprintf('CIRC2: Calculated = %f   Theoretical = %f',...
    tot_band_flux_CIRC2,tot_flux_CIRC2)
%CIRC4 Results
sprintf('CIRC4: Calculated = %f   Theoretical = %f',...
    tot_band_flux_CIRC4,tot_flux_CIRC4)

%Pre-Project Task #2: Change in flux for each band per increase in Kelvin

%Fluxes In Each Band For Original Surface Temp
flux_per_band_old_CIRC2 = band_flux_CIRC2*pi;
flux_per_band_old_CIRC4 = band_flux_CIRC4*pi; 

%Find New Fluxes
syms B_CIRC2(v);syms B_CIRC4(v) 
T_CIRC2 = 307.74+1; %Temperature in Kelvin
T_CIRC4 = 266.55+1; %Temperature in Kelvin
B_CIRC2(v) = 2*h*(c.^2)*(v.^3)*(1/(exp(h*c*v/(k*T_CIRC2))-1));
B_CIRC4(v) = 2*h*(c.^2)*(v.^3)*(1/(exp(h*c*v/(k*T_CIRC4))-1));
B_CIRC2 = matlabFunction(B_CIRC2);
B_CIRC4 = matlabFunction(B_CIRC4);

%CIRC2 Profile
for i = 1:length(bands)-1
    band_flux_CIRC2_new(i) = integral(B_CIRC2, bands(i), bands(i+1));
end

%CIRC4 Profile
for i = 1:length(bands)-1
    band_flux_CIRC4_new(i) = integral(B_CIRC4, bands(i), bands(i+1));
end

%Fluxes In Each Band For New Surface Temp
flux_per_band_new_CIRC2 = band_flux_CIRC2_new*pi;
flux_per_band_new_CIRC4 = band_flux_CIRC4_new*pi;

%Calculate Differences For each Atmospheric Profile
diff_CIRC2 = flux_per_band_new_CIRC2 - flux_per_band_old_CIRC2;
diff_CIRC4 = flux_per_band_new_CIRC4 - flux_per_band_old_CIRC4;


