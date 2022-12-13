%Benjamin Henry 
%AOS527: Atmospheric Radiative Transfer
%Course Project (Part 3: Project Task 1)
%
%
%
%This script handles calculating the TOA and SFC fluxes for an atmopshere
%with each of the 5 gases individually (CO2,H2O,CH4,N2O,O3) and for all 5
%gases together. 

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

%Gas SFC Flux Calculations (Task 1)
%The design for the for loops for all surface flux calculations is that it
%stores the flux per layer of each band in 14 separate columns. The
%resulting arrays progress upwards in layers until Layer N is reached, and
%the final contribution to the arrays is that of the first layer. The for
%loops iterate for each band through all the layers before progressing to
%the next band. The arrays are summed by column to find the total flux per
%band and then completely summed to find total flux per gas.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO2 Calculations
%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.CO2.SFC(:,1))-2
        Temp = data.CIRC2.CO2.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.CO2.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.CO2.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC2.CO2.SFC(j+2,i+2); %transmittance for each layer
        CO2_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.CO2.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            CO2_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.CO2.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_CO2_SFC_band = sum(CO2_flux_layer);
%Total Flux 
tot_CIRC2_CO2_SFC = sum(CO2_flux_layer, 'all');



%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.CO2.SFC(:,1))-2
        Temp = data.CIRC4.CO2.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.CO2.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.CO2.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC4.CO2.SFC(j+2,i+2); %transmittance for each layer
        CO2_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.CO2.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            CO2_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.CO2.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_CO2_SFC_band = sum(CO2_flux_layer);
%Total Flux 
tot_CIRC4_CO2_SFC = sum(CO2_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%H2O Calculations (H2OCTM)
%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.H2OCTM.SFC(:,1))-2
        Temp = data.CIRC2.H2OCTM.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.H2OCTM.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.H2OCTM.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC2.H2OCTM.SFC(j+2,i+2); %transmittance for each layer
        H2O_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.H2OCTM.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            H2O_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.H2OCTM.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_H2OCTM_SFC_band = sum(H2O_flux_layer);
%Total Flux 
tot_CIRC2_H2OCTM_SFC = sum(H2O_flux_layer, 'all');

%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.H2OCTM.SFC(:,1))-2
        Temp = data.CIRC4.H2OCTM.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.H2OCTM.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.H2OCTM.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC4.H2OCTM.SFC(j+2,i+2); %transmittance for each layer
        H2O_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.H2OCTM.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            H2O_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.H2OCTM.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_H2OCTM_SFC_band = sum(H2O_flux_layer);
%Total Flux 
tot_CIRC4_H2OCTM_SFC = sum(H2O_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CH4 Calculations
%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.CH4.SFC(:,1))-2
        Temp = data.CIRC2.CH4.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.CH4.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.CH4.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC2.CH4.SFC(j+2,i+2); %transmittance for each layer
        CH4_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.CH4.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            CH4_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.CH4.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_CH4_SFC_band = sum(CH4_flux_layer);
%Total Flux 
tot_CIRC2_CH4_SFC = sum(CH4_flux_layer, 'all');

%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.CH4.SFC(:,1))-2
        Temp = data.CIRC4.CH4.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.CH4.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.CH4.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC4.CH4.SFC(j+2,i+2); %transmittance for each layer
        CH4_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.CH4.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            CH4_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.CH4.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_CH4_SFC_band = sum(CH4_flux_layer);
%Total Flux 
tot_CIRC4_CH4_SFC = sum(CH4_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%N2O Calculations
%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.N2O.SFC(:,1))-2
        Temp = data.CIRC2.N2O.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.N2O.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.N2O.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC2.N2O.SFC(j+2,i+2); %transmittance for each layer
        N2O_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.N2O.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            N2O_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.N2O.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_N2O_SFC_band = sum(N2O_flux_layer);
%Total Flux 
tot_CIRC2_N2O_SFC = sum(N2O_flux_layer, 'all');

%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.N2O.SFC(:,1))-2
        Temp = data.CIRC4.N2O.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.N2O.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.N2O.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC4.N2O.SFC(j+2,i+2); %transmittance for each layer
        N2O_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.N2O.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            N2O_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.N2O.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_N2O_SFC_band = sum(N2O_flux_layer);
%Total Flux 
tot_CIRC4_N2O_SFC = sum(N2O_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%O3 Calculations
%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.O3.SFC(:,1))-2
        Temp = data.CIRC2.O3.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.O3.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.O3.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC2.O3.SFC(j+2,i+2); %transmittance for each layer
        O3_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.O3.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            O3_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.O3.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_O3_SFC_band = sum(O3_flux_layer);
%Total Flux 
tot_CIRC2_O3_SFC = sum(O3_flux_layer, 'all');

%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.O3.SFC(:,1))-2
        Temp = data.CIRC4.O3.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.O3.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.O3.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC4.O3.SFC(j+2,i+2); %transmittance for each layer
        O3_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.O3.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            O3_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.O3.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_O3_SFC_band = sum(O3_flux_layer);
%Total Flux 
tot_CIRC4_O3_SFC = sum(O3_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%All Gas Calculations
%CIRC2-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC2.ALLGAS.SFC(:,1))-2
        Temp = data.CIRC2.ALLGAS.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.ALLGAS.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC2.ALLGAS.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC2.ALLGAS.SFC(j+2,i+2); %transmittance for each layer
        ALLGAS_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC2.ALLGAS.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            ALLGAS_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC2.ALLGAS.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_ALLGAS_SFC_band = sum(ALLGAS_flux_layer);
%Total Flux 
tot_CIRC2_ALLGAS_SFC = sum(ALLGAS_flux_layer, 'all');

%CIRC4-SFC
for i = 1:length(bands)-1
    for j=1:length(data.CIRC4.ALLGAS.SFC(:,1))-2
        Temp = data.CIRC4.ALLGAS.SFC(j+1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.ALLGAS.SFC(1,1)); %planck fn @ layer 1
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        tb_i = data.CIRC4.ALLGAS.SFC(j+1,i+2); %transmittance for each layer
        tb_iplus1 = data.CIRC4.ALLGAS.SFC(j+2,i+2); %transmittance for each layer
        ALLGAS_flux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iplus1); %Flux Per Layer via Lecture 3 Slides
        if j == length(data.CIRC4.ALLGAS.SFC(:,1))-2
            %Add Layer 1 contribution at the end like in lecture slides
            ALLGAS_flux_layer(j+1,i) = pi*integral(P_layer1,bands(i),...
                bands(i+1))*(1-data.CIRC4.ALLGAS.SFC(1,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_ALLGAS_SFC_band = sum(ALLGAS_flux_layer);
%Total Flux 
tot_CIRC4_ALLGAS_SFC = sum(ALLGAS_flux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gas TOA Flux Calculations (Task 1)
%The design for the for loops for all TOA flux calculations is that it
%stores the flux per layer of each band in 14 separate columns. The
%resulting arrays progress downwards in layers from Layer N-1 until Layer 
%1 is reached, and the final contribution to the arrays is that of the 
%first layer and Layer N as specified in the Lecture 3 notes. The for
%loops iterate for each band through all the layers before progressing to 
%the next band. The arrays are summed by column to find the total flux per
%band and then completely summed to find total flux per gas.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO2 Calculations
%CIRC2-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC2.CO2.TOA(:,1))-1:-1:2
        Temp = data.CIRC2.CO2.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.CO2.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC2.CO2.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC2.CO2.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC2.CO2.TOA(j-1,i+2); %transmittance for each layer
        CO2_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             CO2_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.CO2.TOA(1,i+2));
             %Layer N
             CO2_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.CO2.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_CO2_TOA_band = sum(CO2_TOAflux_layer);
%Total Flux 
tot_CIRC2_CO2_TOA = sum(CO2_TOAflux_layer, 'all');

%CIRC4-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC4.CO2.TOA(:,1))-1:-1:2
        Temp = data.CIRC4.CO2.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.CO2.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC4.CO2.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC4.CO2.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC4.CO2.TOA(j-1,i+2); %transmittance for each layer
        CO2_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             CO2_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.CO2.TOA(1,i+2));
             %Layer N
             CO2_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.CO2.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_CO2_TOA_band = sum(CO2_TOAflux_layer);
%Total Flux 
tot_CIRC4_CO2_TOA = sum(CO2_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%H2OCTM Calculations (H2OCTMCTM)
%CIRC2-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC2.H2OCTM.TOA(:,1))-1:-1:2
        Temp = data.CIRC2.H2OCTM.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.H2OCTM.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC2.H2OCTM.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC2.H2OCTM.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC2.H2OCTM.TOA(j-1,i+2); %transmittance for each layer
        H2OCTM_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             H2OCTM_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.H2OCTM.TOA(1,i+2));
             %Layer N
             H2OCTM_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.H2OCTM.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_H2OCTM_TOA_band = sum(H2OCTM_TOAflux_layer);
%Total Flux 
tot_CIRC2_H2OCTM_TOA = sum(H2OCTM_TOAflux_layer, 'all');

%CIRC4-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC4.H2OCTM.TOA(:,1))-1:-1:2
        Temp = data.CIRC4.H2OCTM.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.H2OCTM.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC4.H2OCTM.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC4.H2OCTM.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC4.H2OCTM.TOA(j-1,i+2); %transmittance for each layer
        H2OCTM_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             H2OCTM_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.H2OCTM.TOA(1,i+2));
             %Layer N
             H2OCTM_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.H2OCTM.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_H2OCTM_TOA_band = sum(H2OCTM_TOAflux_layer);
%Total Flux 
tot_CIRC4_H2OCTM_TOA = sum(H2OCTM_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CH4 Calculations
%CIRC2-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC2.CH4.TOA(:,1))-1:-1:2
        Temp = data.CIRC2.CH4.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.CH4.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC2.CH4.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC2.CH4.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC2.CH4.TOA(j-1,i+2); %transmittance for each layer
        CH4_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             CH4_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.CH4.TOA(1,i+2));
             %Layer N
             CH4_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.CH4.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_CH4_TOA_band = sum(CH4_TOAflux_layer);
%Total Flux 
tot_CIRC2_CH4_TOA = sum(CH4_TOAflux_layer, 'all');

%CIRC4-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC4.CH4.TOA(:,1))-1:-1:2
        Temp = data.CIRC4.CH4.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.CH4.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC4.CH4.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC4.CH4.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC4.CH4.TOA(j-1,i+2); %transmittance for each layer
        CH4_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             CH4_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.CH4.TOA(1,i+2));
             %Layer N
             CH4_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.CH4.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_CH4_TOA_band = sum(CH4_TOAflux_layer);
%Total Flux 
tot_CIRC4_CH4_TOA = sum(CH4_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%N2O Calculations
%CIRC2-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC2.N2O.TOA(:,1))-1:-1:2
        Temp = data.CIRC2.N2O.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.N2O.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC2.N2O.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC2.N2O.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC2.N2O.TOA(j-1,i+2); %transmittance for each layer
        N2O_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             N2O_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.N2O.TOA(1,i+2));
             %Layer N
             N2O_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.N2O.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_N2O_TOA_band = sum(N2O_TOAflux_layer);
%Total Flux 
tot_CIRC2_N2O_TOA = sum(N2O_TOAflux_layer, 'all');

%CIRC4-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC4.N2O.TOA(:,1))-1:-1:2
        Temp = data.CIRC4.N2O.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.N2O.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC4.N2O.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC4.N2O.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC4.N2O.TOA(j-1,i+2); %transmittance for each layer
        N2O_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             N2O_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.N2O.TOA(1,i+2));
             %Layer N
             N2O_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.N2O.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_N2O_TOA_band = sum(N2O_TOAflux_layer);
%Total Flux 
tot_CIRC4_N2O_TOA = sum(N2O_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%O3 Calculations
%CIRC2-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC2.O3.TOA(:,1))-1:-1:2
        Temp = data.CIRC2.O3.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC2.O3.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC2.O3.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC2.O3.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC2.O3.TOA(j-1,i+2); %transmittance for each layer
        O3_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             O3_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.O3.TOA(1,i+2));
             %Layer N
             O3_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.O3.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_O3_TOA_band = sum(O3_TOAflux_layer);
%Total Flux 
tot_CIRC2_O3_TOA = sum(O3_TOAflux_layer, 'all');

%CIRC4-TOA
for i = 1:length(bands)-1
    for j=length(data.CIRC4.O3.TOA(:,1))-1:-1:2
        Temp = data.CIRC4.O3.TOA(j-1,1); %temperature of each layer
        P_layer = B(v,Temp); %planck fn for each layer
        P_layer1 = B(v,data.CIRC4.O3.TOA(1,1)); %planck fn @ layer 1
        P_layerN = B(v,data.CIRC4.O3.TOA(end,1)); %planck fn @ layer N
        P_layer = matlabFunction(P_layer);
        P_layer1 = matlabFunction(P_layer1);
        P_layerN = matlabFunction(P_layerN);
        tb_i = data.CIRC4.O3.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC4.O3.TOA(j-1,i+2); %transmittance for each layer
        O3_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             O3_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.O3.TOA(1,i+2));
             %Layer N
             O3_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.O3.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_O3_TOA_band = sum(O3_TOAflux_layer);
%Total Flux 
tot_CIRC4_O3_TOA = sum(O3_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ALLGAS Calculations

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
        tb_i = data.CIRC2.ALLGAS.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC2.ALLGAS.TOA(j-1,i+2); %transmittance for each layer
        ALLGAS_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             ALLGAS_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC2.ALLGAS.TOA(1,i+2));
             %Layer N
             ALLGAS_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC2.ALLGAS.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC2_ALLGAS_TOA_band = sum(ALLGAS_TOAflux_layer);
%Total Flux 
tot_CIRC2_ALLGAS_TOA = sum(ALLGAS_TOAflux_layer, 'all');

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
        tb_i = data.CIRC4.ALLGAS.TOA(j,i+2); %transmittance for each layer
        tb_iminus1 = data.CIRC4.ALLGAS.TOA(j-1,i+2); %transmittance for each layer
        ALLGAS_TOAflux_layer(j,i) = pi*integral(P_layer,bands(i),bands(i+1))*...
            (tb_i-tb_iminus1);%Flux Per Layer via Lecture 3 Slides
        if j == 2
            %Add Layer 1 & N contributions at the end 
            %Layer 1
             ALLGAS_TOAflux_layer(1,i) = pi*integral(P_layer1,bands(i), ...
                 bands(i+1))*(data.CIRC4.ALLGAS.TOA(1,i+2));
             %Layer N
             ALLGAS_TOAflux_layer(54,i) = ...
                 pi*integral(P_layerN,bands(i), ...
                 bands(i+1))*(1-data.CIRC4.ALLGAS.TOA(54,i+2));
        end
    end
end
%Total Flux per band
tot_CIRC4_ALLGAS_TOA_band = sum(ALLGAS_TOAflux_layer);
%Total Flux 
tot_CIRC4_ALLGAS_TOA = sum(ALLGAS_TOAflux_layer, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

