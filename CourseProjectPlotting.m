%Benjamin Henry
%AOS527:Atmospheric Radiative Transfer
%Course Project (Plotting)
%
%
%
%This script handles plotting all of the results. This script is separate
%from the main computations in order to cut down run times (run time of
%main scripts, CourseProjectPt3 & CourseProjectPt4, is near 15-20 minutes).
%It is important to note that the results of the main script are used in
%order to plot here, so they must be run BEFORE plotting here. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TASK 1 PLOTTING: SFC & TOA FLUXES

x = 1:1:length(bands)-1; %x-coordinate for all plots

%CIRC2 SFC FLUXES
figure(1)
subplot(2,2,1)
plot(x,tot_CIRC2_CO2_SFC_band,'.-')
hold on 
plot(x,tot_CIRC2_H2OCTM_SFC_band,'.-')
plot(x,tot_CIRC2_CH4_SFC_band,'.-')
plot(x,tot_CIRC2_N2O_SFC_band,'.-')
plot(x,tot_CIRC2_O3_SFC_band,'.-')
plot(x,tot_CIRC2_ALLGAS_SFC_band,'.-')
title('CIRC2 SFC Fluxes for Each Band','.-')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('CO2 (Total = 98.7 W/m^2)','H2OCTM (Total = 457.8 W/m^2)',...
    'CH4 (Total = 7.2 W/m^2)','N2O (Total = 8.4 W/m^2)',...
    'O3 (Total = 7.2 W/m^2)','ALLGAS (Total = 490.3 W/m^2)')
hold off


%CIRC2 TOA FLUXES
subplot(2,2,2)
plot(x,tot_CIRC2_CO2_TOA_band,'.-')
hold on 
plot(x,tot_CIRC2_H2OCTM_TOA_band,'.-')
plot(x,tot_CIRC2_CH4_TOA_band,'.-')
plot(x,tot_CIRC2_N2O_TOA_band,'.-')
plot(x,tot_CIRC2_O3_TOA_band,'.-')
plot(x,tot_CIRC2_ALLGAS_TOA_band,'.-')
title('CIRC2 TOA Fluxes for Each Band','.-')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('CO2 (Total = 452.2 W/m^2)','H2OCTM (Total = 335.3 W/m^2)',...
    'CH4 (Total = 499.9 W/m^2)','N2O (Total = 499.3 W/m^2)',...
    'O3 (Total = 492.9 W/m^2)','ALLGAS (Total = 291.0 W/m^2)')
hold off

%CIRC4 SFC FLUXES
subplot(2,2,3)
plot(x,tot_CIRC4_CO2_SFC_band,'.-')
hold on 
plot(x,tot_CIRC4_H2OCTM_SFC_band,'.-')
plot(x,tot_CIRC4_CH4_SFC_band,'.-')
plot(x,tot_CIRC4_N2O_SFC_band,'.-')
plot(x,tot_CIRC4_O3_SFC_band,'.-')
plot(x,tot_CIRC4_ALLGAS_SFC_band,'.-')
title('CIRC4 SFC Fluxes for Each Band','.-')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('CO2 (Total = 48.6 W/m^2)','H2OCTM (Total = 135.6 W/m^2)',...
    'CH4 (Total = 3.6 W/m^2)','N2O (Total = 4.8 W/m^2)',...
    'O3 (Total = 5.3 W/m^2)','ALLGAS (Total = 176.0 W/m^2)')
hold off

%CIRC4 TOA FLUXES
subplot(2,2,4)
plot(x,tot_CIRC4_CO2_TOA_band,'.-')
hold on 
plot(x,tot_CIRC4_H2OCTM_TOA_band,'.-')
plot(x,tot_CIRC4_CH4_TOA_band,'.-')
plot(x,tot_CIRC4_N2O_TOA_band,'.-')
plot(x,tot_CIRC4_O3_TOA_band,'.-')
plot(x,tot_CIRC4_ALLGAS_TOA_band,'.-')
title('CIRC4 TOA Fluxes for Each Band','.-')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('CO2 (Total = 264.3 W/m^2)','H2OCTM (Total = 251.8 W/m^2)',...
    'CH4 (Total = 281.7 W/m^2)','N2O (Total = 281.6 W/m^2)',...
    'O3 (Total = 277.4 W/m^2)','ALLGAS (Total = 226.3 W/m^2)')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TASK 2 PLOTTING: WATER VAPOR CONTINUUM

%H2OCTM ONLY
figure(2)
plot(x,tot_CIRC2_H2OCTM_SFC_band,'.-')
hold on
plot(x,tot_CIRC2_H2OCTM_TOA_band,'.-')
plot(x,tot_CIRC4_H2OCTM_SFC_band,'.-')
plot(x,tot_CIRC4_H2OCTM_TOA_band,'.-')
title('H2OCTM Fluxes')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('CIRC2 SFC','CIRC2 TOA','CIRC4 SFC','CIRC4 TOA')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TASK 3 PLOTTING: CO22x Analysis

%CIRC2
figure(3)
subplot(2,1,1)
plot(x,tot_CIRC2_CO22x_SFC_band,'.-')
hold on
plot(x,tot_CIRC2_CO22x_TOA_band,'.-')
plot(x,tot_CIRC2_CO2_SFC_band,'.-')
plot(x,tot_CIRC2_CO2_TOA_band,'.-')
title('CIRC2 Radiative Fluxes For CO2 vs. CO22x')
ylabel('Flux (W/m^2')
xlabel('Band #')
legend('SFC CO22x (Total = 113.8 W/m^2)','TOA CO22x (Total = 444.6 W/m^2)'...
    ,'SFC CO2 (Total = 98.7 W/m^2)', 'TOA CO2 (Total = 452.2 W/m^2)')
hold off

%CIRC4
subplot(2,1,2)
plot(x,tot_CIRC4_CO22x_SFC_band,'.-')
hold on
plot(x,tot_CIRC4_CO22x_TOA_band,'.-')
plot(x,tot_CIRC4_CO2_SFC_band,'.-')
plot(x,tot_CIRC4_CO2_TOA_band,'.-')
title('CIRC4 Radiative Fluxes For CO2 vs. CO22x')
ylabel('Flux (W/m^2')
xlabel('Band #')
legend('SFC CO22x (Total = 53.6 W/m^2)','TOA CO22x (Total = 262.8 W/m^2)'...
    ,'SFC CO2 (Total = 48.6 W/m^2)', 'TOA CO2 (Total = 264.3 W/m^2)')
hold off

%CIRC2 ALLGAS Without and With Water Vapor Continuum
figure(4)
plot(x, tot_CIRC2_ALLGAS_SFC_band_NOCTM)
hold on
plot(x, tot_CIRC2_ALLGAS_SFC_bandCTM)
plot(x, tot_CIRC2_ALLGAS_TOA_band_NOCTM)
plot(x, tot_CIRC2_ALLGAS_TOA_bandCTM)
title('CIRC2 Radiative Fluxes For All Gases with/without H2O Continuum under 2xCO2')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('SFC no CTM (Total = 441.4 W/m^2)','SFC w/ CTM (Total = 494.0 W/m^2)'...
    ,'TOA no CTM (Total = 262.3 W/m^2)', 'TOA w/ CTM (Total = 215.5 W/m^2)')
hold off

%CIRC4 ALLGAS Without and With Water Vapor Continuum
figure(5)
plot(x, tot_CIRC4_ALLGAS_SFC_band_NOCTM)
hold on
plot(x, tot_CIRC4_ALLGAS_SFC_bandCTM)
plot(x, tot_CIRC4_ALLGAS_TOA_band_NOCTM)
plot(x, tot_CIRC4_ALLGAS_TOA_bandCTM)
title('CIRC4 Radiative Fluxes For All Gases with/without H2O Continuum under 2xCO2')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('SFC no CTM (Total = 172.9 W/m^2)','SFC w/ CTM (Total = 178.4 W/m^2)'...
    ,'TOA no CTM (Total = 211.5 W/m^2)', 'TOA w/ CTM (Total = 209.1 W/m^2)')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TASK 4 PLOTTING: Global Warming Temperature Profiles

%Plotting the Changes in TOA FLUXES
figure(6)
plot(x,CIRC2_ALLGAS_TOA_CHANGE_band)
hold on 
plot(x,CIRC4_ALLGAS_TOA_CHANGE_band)
plot(x, diff_CIRC2)
plot(x, diff_CIRC4)
title('Difference in TOA Fluxes')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('CIRC2 (W/ ATM)','CIRC4 (W/ ATM)','CIRC2 (NO ATM)', 'CIRC4 (NO ATM)')
hold off

%Plotting the New TOA FLUXES
figure(7)
plot(x,newtot_CIRC2_ALLGAS_TOA_band)
hold on 
plot(x,newtot_CIRC4_ALLGAS_TOA_band)
plot(x,flux_per_band_new_CIRC2)
plot(x,flux_per_band_new_CIRC4)
title('TOA Fluxes from New Temperature Profile')
ylabel('Flux (W/m^2)')
xlabel('Band #')
legend('CIRC2 (W/ ATM)(Total = 294.6 W/m^2)',...
    'CIRC4 (W/ ATM)(Total = 229.0 W/m^2)',...
    'CIRC2 (NO ATM)(Total = 512.4 W/m^2)',...
    'CIRC4 (NO ATM)(Total = 289.1 W/m^2)')
hold off
