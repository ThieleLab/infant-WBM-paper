function InputDataMetabolitesBC = integrateNBSMetaboliteConcentrations(NB_data)
% This function converts newborn screening metabolite concentration into a
% input format such that they can be integrated using the
% physiologicalConstraintsHMDBbased function. 
% metabolite concentrations have to be given in uM

% INPUT
% NB_data                   table row with metabolite concentration (umol/L)
%                           starting in column 6
%
% OUTPUT
% InputDataMetabolitesBC    cell with min and max values for metabolite
%                           concentrations measured in (umol/L)
%
% Elaine Zaunseder, 2022


%Load default parameters
NB_data_array=table2array(NB_data);

%Load variation coefficients given by Heidelberg NBS center
%This should be updated center-specific for other NBS centers
variation_data = readtable('variation_coefficient.xlsx');

%iterate over input metabolites 
for j= 6: (length(NB_data_array)-1) %starts at the first metabolite concentration
   %load metabolite concentrations
   metabolite=cell2mat(NB_data.Properties.VariableNames(j));
   concentration= NB_data_array(1,j); %micromol(umol)/Liter
   
   %find index of variation coefficient
   ind=find(strcmp(variation_data.Properties.VariableNames,NB_data.Properties.VariableNames(j)));
   variation = table2array(variation_data(1,ind));
   
   % for some metabolites the std of the internal standard is very large,
   % therefore we limit the variation to 10% in these cases
   if variation > 0.1
       variation=0.1;
   end

   % calculate flux boundaries - EZ 24.2.2023
   fluxMin = concentration *(1-variation); %in micromol(umol)/L
   fluxMax = concentration *(1+variation); %in micromol(umol)/L
   

  if fluxMin < 0
      fluxMin = 0;
  end
   
   %formatting for output file
   fluxMin = cellstr(num2str(fluxMin));
   fluxMin = regexprep(fluxMin,' ','');
   for i = 1 : size(fluxMin,1)
        fluxMin{i,1} = (fluxMin(i,1));
   end
   fluxMax = cellstr(num2str(fluxMax));     
   fluxMax = regexprep(fluxMax,' ','');
   for i = 1 : size(fluxMax,1)
        fluxMax{i,1} = (fluxMax(i,1));
   end
   
   %Update Input data table 
   InputDataMetabolitesBC(j-5,:) = [metabolite   fluxMin  fluxMax];
   
end