function results = makeActivityModel(baseModel, amplitude)
% makeActivityModel
% Loads the activity values (values between 0 and 1) from the file. And
% multiply with amplitude (kcal/kg) to get the energy expenditure for
% activity at each timepoint.
%
%   baseModel              a list of values between 0 and 1. first column
%                          contains the days starting from 1
%   amplitude              a value in kcal/kg to scale the model with
%   results                a list of kcal/kg, one value for each day
%
%   Avlant Nilsson, 2016-05-16
%   This published acitivity model and code stems from Nilsson et al
%   (https://doi.org/10.1038/s41540-017-0004-5)

 %Exchange reactions, Metabolites
    results = load(baseModel);
    results(:,2) = results(:,2) * amplitude;
    results(:,1) = [];
end

