% Benjamin Henry 
% AOS527: Atmospheric Radiative Transfer 
% Course Project(Part 1: Data Collection)
%
%
% The script below reads in all of the necessary data files and saves them
% locally in a structure so that the data files do not need to be read
% every time the calculation scripts are run, cutting down run times. 

%Adding the Proper Directory to Retreive 
%addpath('/Users/bshenry/Documents/MATLAB/AOS527/CourseProject')

%Creating Matrices for differing gasses, profiles, and TOA vs. SFC.
comp = {'2xCO2', 'ALLGAS', 'CH4', 'CO2', 'H2O', 'H2OCTM', 'N2O', 'O3'};
prof = {'CIRC1', 'CIRC2', 'CIRC4'};
sfc_or_toa = {'SFC','TOA'};

%Structure To Store the Matrices
data = struct();

%CIRC2 & CIRC4 Data 
for i=1:length(comp) %for all gasses
    for j = 1:2 %for either the sfc or toa
        for l = 2:3 %for CIRC2 and CIRC4 (not CIRC1)
            txtfile = char(strcat(prof(l),'_LAYER_',...
                comp(i),'_',sfc_or_toa(j),'.txt'));
            if strcmp(comp(i),'2xCO2')==1
                 % 2xCO2 cannot be added as field name bc it starts with number
                data.(prof{l}).CO22x.(sfc_or_toa{j}) = readmatrix(txtfile);
            else
                % add matrix into the structure
                data.(prof{l}).(comp{i}).(sfc_or_toa{j}) = readmatrix(txtfile); 
        
            end
        end
    end
end
% getting CIRC1 data
data.(prof{1}).(comp{5}).(sfc_or_toa{1}) = readmatrix('CIRC1_LAYER_H2O_SFC.txt');
data.(prof{1}).(comp{5}).(sfc_or_toa{2}) = readmatrix('CIRC1_LAYER_H2O_TOA.txt');
