%Benjamin Henry 
%AOS527: Atmospheric Radiative Transfer
%Course Project (Part 3: Project Task 3)
%
%
%
%This script handles calculating the TOA and SFC fluxes for an atmopshere
%double the amount of CO22x with no other gases included, and with the 4
%other gases included.

%Preliminaries
bands = [0 160 560 630 700 800 900 990 1070 1200 1400 ...
    2200 2250 2400 3000];  %Bands in Wavenumber limits (cm^-1)
bands = bands*100; %Conversion from cm^-1 to m^-1
h = 6.63*10.^(-34); %Plank Constant (J*s)
c = 3*10.^8; %Speed of Light (m/s)
k = 1.38*10.^(-23); %Boltzmann Constant (m^2*kg/s^2*K)
syms B(v,T)
B(v,T) = 2*h*(c.^2)*(v.^3)*(1/(exp(h*c*v/(k*T))-1));
B = matlabFunction(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO22x Doubling Calculations (Task 3)
%Identical Design to that found in CourseProjectPt3, in that the flux for
%each band is found by taking each band, iterating through all of the
%layers, and then moving on to the next band. The total flux per band is
%found, followed by the total overall flux.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO22x Only Calculations

%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.CO22x.SFC(:,1))-2
        Temp = data.CIRC2.CO22x.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.CO22x.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.CO22x.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC2.CO22x.SFC(j+2,i+2); %transmittance for each layer
        CO22x_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.CO22x.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            CO22x_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.CO22x.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_CO22x_SFC_band = sum(CO22x_flux_layer);
%Total Flux 
tot_CIRC2_CO22x_SFC = sum(CO22x_flux_layer, 'all');

%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.CO22x.SFC(:,1))-2
        Temp = data.CIRC4.CO22x.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.CO22x.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.CO22x.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC4.CO22x.SFC(j+2,i+2); %transmittance for each layer
        CO22x_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.CO22x.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            CO22x_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.CO22x.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_CO22x_SFC_band = sum(CO22x_flux_layer);
%Total Flux 
tot_CIRC4_CO22x_SFC = sum(CO22x_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CIRC2-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC2.CO22x.TOA(:,1))-1:-1:2
        Temp = data.CIRC2.CO22x.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.CO22x.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC2.CO22x.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC2.CO22x.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC2.CO22x.TOA(j-1,i+2); %transmittance for each layer
        CO22x_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             CO22x_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.CO22x.TOA(1,i+2));
             %Layer N
             CO22x_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.CO22x.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_CO22x_TOA_band = sum(CO22x_TOAflux_layer);
%Total Flux 
tot_CIRC2_CO22x_TOA = sum(CO22x_TOAflux_layer, 'all');

%CIRC4-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC4.CO22x.TOA(:,1))-1:-1:2
        Temp = data.CIRC4.CO22x.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.CO22x.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC4.CO22x.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC4.CO22x.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC4.CO22x.TOA(j-1,i+2); %transmittance for each layer
        CO22x_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             CO22x_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.CO22x.TOA(1,i+2));
             %Layer N
             CO22x_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.CO22x.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_CO22x_TOA_band = sum(CO22x_TOAflux_layer);
%Total Flux 
tot_CIRC4_CO22x_TOA = sum(CO22x_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO22x Doubling Calculations (Task 3)
%Identical Design to that found in CourseProjectPt3, in that the flux for
%each band is found by taking each band, iterating through all of the
%layers, and then moving on to the next band. The total flux per band is
%found, followed by the total overall flux.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO22x All Gas SFC Calculations (With Water Vapor Continuum)

%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.CO22x.SFC(:,1))-2
        Temp = data.CIRC2.CO22x.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.CO2.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.CO22x.SFC(j+1,i+2)*data.CIRC2.CH4.SFC(j+1,i+2)...
            *data.CIRC2.H2OCTM.SFC(j+1,i+2)*data.CIRC2.N2O.SFC(j+1,i+2)...
            *data.CIRC2.O3.SFC(j+1,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        tb_iplus1 = data.CIRC2.CO22x.SFC(j+2,i+2)*data.CIRC2.CH4.SFC(j+2,...
            i+2)*data.CIRC2.H2OCTM.SFC(j+2,i+2)*data.CIRC2.N2O.SFC(j+2,i+2)...
            *data.CIRC2.O3.SFC(j+2,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        ALLGAS_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.CO2.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            ALLGAS_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.CO22x.SFC(1,i+2)*data.CIRC2.CH4.SFC(1,...
            i+2)*data.CIRC2.H2OCTM.SFC(1,i+2)*data.CIRC2.N2O.SFC(1,i+2)...
            *data.CIRC2.O3.SFC(1,i+2));
            
        end
    end
end
%Total Flux per band
tot_CIRC2_ALLGAS_SFC_bandCTM = sum(ALLGAS_flux_layer);
%Total Flux 
tot_CIRC2_ALLGAS_SFCCTM = sum(ALLGAS_flux_layer, 'all');

%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.CO22x.SFC(:,1))-2
        Temp = data.CIRC4.CO22x.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.CO2.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.CO22x.SFC(j+1,i+2)*data.CIRC4.CH4.SFC(j+1,i+2)...
            *data.CIRC4.H2OCTM.SFC(j+1,i+2)*data.CIRC4.N2O.SFC(j+1,i+2)...
            *data.CIRC4.O3.SFC(j+1,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        tb_iplus1 = data.CIRC4.CO22x.SFC(j+2,i+2)*data.CIRC4.CH4.SFC(j+2,...
            i+2)*data.CIRC4.H2OCTM.SFC(j+2,i+2)*data.CIRC4.N2O.SFC(j+2,i+2)...
            *data.CIRC4.O3.SFC(j+2,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        ALLGAS_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.CO2.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            ALLGAS_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.CO22x.SFC(1,i+2)*data.CIRC4.CH4.SFC(1,...
            i+2)*data.CIRC4.H2OCTM.SFC(1,i+2)*data.CIRC4.N2O.SFC(1,i+2)...
            *data.CIRC4.O3.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_ALLGAS_SFC_bandCTM = sum(ALLGAS_flux_layer);
%Total Flux 
tot_CIRC4_ALLGAS_SFCCTM = sum(ALLGAS_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO22x All Gas SFC Calculations (Without Water Vapor Continuum)

%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.CO22x.SFC(:,1))-2
        Temp = data.CIRC2.CO22x.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.CO2.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.CO22x.SFC(j+1,i+2)*data.CIRC2.CH4.SFC(j+1,i+2)...
            *data.CIRC2.H2O.SFC(j+1,i+2)*data.CIRC2.N2O.SFC(j+1,i+2)...
            *data.CIRC2.O3.SFC(j+1,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H2O, CH4, N20, 03)
        tb_iplus1 = data.CIRC2.CO22x.SFC(j+2,i+2)*data.CIRC2.CH4.SFC(j+2,...
            i+2)*data.CIRC2.H2O.SFC(j+2,i+2)*data.CIRC2.N2O.SFC(j+2,i+2)...
            *data.CIRC2.O3.SFC(j+2,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H2O, CH4, N20, 03)
        ALLGAS_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.CO2.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            ALLGAS_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.CO22x.SFC(1,i+2)*data.CIRC2.CH4.SFC(1,...
            i+2)*data.CIRC2.H2O.SFC(1,i+2)*data.CIRC2.N2O.SFC(1,i+2)...
            *data.CIRC2.O3.SFC(1,i+2));
            
        end
    end
end
%Total Flux per band
tot_CIRC2_ALLGAS_SFC_band_NOCTM = sum(ALLGAS_flux_layer);
%Total Flux 
tot_CIRC2_ALLGAS_SFC_NOCTM = sum(ALLGAS_flux_layer, 'all');

%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.CO22x.SFC(:,1))-2
        Temp = data.CIRC4.CO22x.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.CO2.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.CO22x.SFC(j+1,i+2)*data.CIRC4.CH4.SFC(j+1,i+2)...
            *data.CIRC4.H2O.SFC(j+1,i+2)*data.CIRC4.N2O.SFC(j+1,i+2)...
            *data.CIRC4.O3.SFC(j+1,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H2O, CH4, N20, 03)
        tb_iplus1 = data.CIRC4.CO22x.SFC(j+2,i+2)*data.CIRC4.CH4.SFC(j+2,...
            i+2)*data.CIRC4.H2O.SFC(j+2,i+2)*data.CIRC4.N2O.SFC(j+2,i+2)...
            *data.CIRC4.O3.SFC(j+2,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H2O, CH4, N20, 03)
        ALLGAS_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.CO2.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            ALLGAS_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.CO22x.SFC(1,i+2)*data.CIRC4.CH4.SFC(1,...
            i+2)*data.CIRC4.H2O.SFC(1,i+2)*data.CIRC4.N2O.SFC(1,i+2)...
            *data.CIRC4.O3.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_ALLGAS_SFC_band_NOCTM = sum(ALLGAS_flux_layer);
%Total Flux 
tot_CIRC4_ALLGAS_SFC_NOCTM = sum(ALLGAS_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO22x All Gas TOA Calculations (With Water Vapor Continuum)

%CIRC2-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC2.ALLGAS.TOA(:,1))-1:-1:2
        Temp = data.CIRC2.ALLGAS.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.ALLGAS.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC2.ALLGAS.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC2.CO22x.SFC(j,i+2)*data.CIRC2.CH4.SFC(j,i+2)...
            *data.CIRC2.H2OCTM.SFC(j,i+2)*data.CIRC2.N2O.SFC(j,i+2)...
            *data.CIRC2.O3.SFC(j,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        tb_iminus1 = data.CIRC2.CO22x.SFC(j-1,i+2)*data.CIRC2.CH4.SFC(j-1,i+2)...
            *data.CIRC2.H2OCTM.SFC(j-1,i+2)*data.CIRC2.N2O.SFC(j-1,i+2)...
            *data.CIRC2.O3.SFC(j-1,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        ALLGAS_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             ALLGAS_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.CO22x.SFC(1,i+2)*data.CIRC2.CH4.SFC(1,i+2)...
            *data.CIRC2.H2OCTM.SFC(1,i+2)*data.CIRC2.N2O.SFC(1,i+2)...
            *data.CIRC2.O3.SFC(1,i+2));
             %Layer N
             ALLGAS_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.CO22x.SFC(54,i+2)*data.CIRC2.CH4.SFC(54,i+2)...
            *data.CIRC2.H2OCTM.SFC(54,i+2)*data.CIRC2.N2O.SFC(54,i+2)...
            *data.CIRC2.O3.SFC(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_ALLGAS_TOA_bandCTM = sum(ALLGAS_TOAflux_layer);
%Total Flux 
tot_CIRC2_ALLGAS_TOACTM = sum(ALLGAS_TOAflux_layer, 'all');

%CIRC4-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC4.ALLGAS.TOA(:,1))-1:-1:2
        Temp = data.CIRC4.ALLGAS.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.ALLGAS.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC4.ALLGAS.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC4.CO22x.SFC(j,i+2)*data.CIRC4.CH4.SFC(j,i+2)...
            *data.CIRC4.H2OCTM.SFC(j,i+2)*data.CIRC4.N2O.SFC(j,i+2)...
            *data.CIRC4.O3.SFC(j,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        tb_iminus1 = data.CIRC4.CO22x.SFC(j-1,i+2)*data.CIRC4.CH4.SFC(j-1,i+2)...
            *data.CIRC4.H2OCTM.SFC(j-1,i+2)*data.CIRC4.N2O.SFC(j-1,i+2)...
            *data.CIRC4.O3.SFC(j-1,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        ALLGAS_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             ALLGAS_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.CO22x.SFC(1,i+2)*data.CIRC4.CH4.SFC(1,i+2)...
            *data.CIRC4.H2OCTM.SFC(1,i+2)*data.CIRC4.N2O.SFC(1,i+2)...
            *data.CIRC4.O3.SFC(1,i+2));
             %Layer N
             ALLGAS_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.CO22x.SFC(54,i+2)*data.CIRC4.CH4.SFC(54,i+2)...
            *data.CIRC4.H2OCTM.SFC(54,i+2)*data.CIRC4.N2O.SFC(54,i+2)...
            *data.CIRC4.O3.SFC(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_ALLGAS_TOA_bandCTM = sum(ALLGAS_TOAflux_layer);
%Total Flux 
tot_CIRC4_ALLGAS_TOACTM = sum(ALLGAS_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO22x All Gas TOA Calculations (Without Water Vapor Continuum)

%CIRC2-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC2.ALLGAS.TOA(:,1))-1:-1:2
        Temp = data.CIRC2.ALLGAS.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.ALLGAS.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC2.ALLGAS.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC2.CO22x.SFC(j,i+2)*data.CIRC2.CH4.SFC(j,i+2)...
            *data.CIRC2.H2O.SFC(j,i+2)*data.CIRC2.N2O.SFC(j,i+2)...
            *data.CIRC2.O3.SFC(j,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        tb_iminus1 = data.CIRC2.CO22x.SFC(j-1,i+2)*data.CIRC2.CH4.SFC(j-1,i+2)...
            *data.CIRC2.H2O.SFC(j-1,i+2)*data.CIRC2.N2O.SFC(j-1,i+2)...
            *data.CIRC2.O3.SFC(j-1,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        ALLGAS_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             ALLGAS_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.CO22x.SFC(1,i+2)*data.CIRC2.CH4.SFC(1,i+2)...
            *data.CIRC2.H2O.SFC(1,i+2)*data.CIRC2.N2O.SFC(1,i+2)...
            *data.CIRC2.O3.SFC(1,i+2));
             %Layer N
             ALLGAS_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.CO22x.SFC(54,i+2)*data.CIRC2.CH4.SFC(54,i+2)...
            *data.CIRC2.H2O.SFC(54,i+2)*data.CIRC2.N2O.SFC(54,i+2)...
            *data.CIRC2.O3.SFC(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_ALLGAS_TOA_band_NOCTM = sum(ALLGAS_TOAflux_layer);
%Total Flux 
tot_CIRC2_ALLGAS_TOA_NOCTM = sum(ALLGAS_TOAflux_layer, 'all');

%CIRC4-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC4.ALLGAS.TOA(:,1))-1:-1:2
        Temp = data.CIRC4.ALLGAS.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.ALLGAS.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC4.ALLGAS.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC4.CO22x.SFC(j,i+2)*data.CIRC4.CH4.SFC(j,i+2)...
            *data.CIRC4.H2O.SFC(j,i+2)*data.CIRC4.N2O.SFC(j,i+2)...
            *data.CIRC4.O3.SFC(j,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        tb_iminus1 = data.CIRC4.CO22x.SFC(j-1,i+2)*data.CIRC4.CH4.SFC(j-1,i+2)...
            *data.CIRC4.H2O.SFC(j-1,i+2)*data.CIRC4.N2O.SFC(j-1,i+2)...
            *data.CIRC4.O3.SFC(j-1,i+2); %transmittance for each layer 
        % found by multiplying the transmittance of each individual gas 
        % together (CO22x, H20CTM, CH4, N20, 03)
        ALLGAS_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             ALLGAS_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.CO22x.SFC(1,i+2)*data.CIRC4.CH4.SFC(1,i+2)...
            *data.CIRC4.H2O.SFC(1,i+2)*data.CIRC4.N2O.SFC(1,i+2)...
            *data.CIRC4.O3.SFC(1,i+2));
             %Layer N
             ALLGAS_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.CO22x.SFC(54,i+2)*data.CIRC4.CH4.SFC(54,i+2)...
            *data.CIRC4.H2O.SFC(54,i+2)*data.CIRC4.N2O.SFC(54,i+2)...
            *data.CIRC4.O3.SFC(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_ALLGAS_TOA_band_NOCTM = sum(ALLGAS_TOAflux_layer);
%Total Flux 
tot_CIRC4_ALLGAS_TOA_NOCTM = sum(ALLGAS_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
