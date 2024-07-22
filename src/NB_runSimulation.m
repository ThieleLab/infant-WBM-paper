function simulationResults = NB_runSimulation(model,simSettings)
% NB_runSimulation(model,simSettings, resultsPATH)
%
% Simulates growth of an infant by solving FBA problems for each day.
% And adapting the weight based on the predicted growth rate from the
% previous day. 

%   model                     The WBM model with which we are starting
%   simSettings               The simulation settings
%   resultsPATH               Path for saving the results
%
%   simulationResults          A struct containing results of the
%                              simulation
%       .timePoints            A list of days (that were simulated).
%       .weight                A list of weights for these days (in gram).
%       .v                     The simulated fluxes
%       .rxns                  The model's reactions names

%   Avlant Nilsson, 2016-05-16
%   Elaine Zaunseder, 2022-2023
%
    modelOri=model;
    
    %get weight table for gender 
    if strcmp(model.sex,'male')
        initialWeight=3.3;
        weights=readtable('/results/3_2_Results/male_growth_models/reference_data/baby_boy_weight.csv');
        T= readtable('dummy_normal_male.xlsx');
    elseif strcmp(model.sex, 'female')
        load('models/babyGirl.mat')
        initialWeight=3.2;
        weights=readtable('/results/3_2_Results/female_growth_models/reference_data/baby_girl_weight.csv');
        T= readtable('dummy_normal_female.xlsx');
    else
        fprintf('Please enter a valid gender, these are "male" or "female"!')
        return
    end
    
    %Set initial weight based on age
    if simSettings.start==1
        weight=initialWeight;
    elseif simSettings.start <=180
        weight=round(weights{simSettings.start-1,1}/1000,4);
    else
        fprintf('Please provide an age in the range 1-180 days!')
        return
    end
    
    %Time steps
    simPoints = getSimPoints(simSettings);
    %Allocate memory
    weights = (zeros(length(simPoints), 1)+weight);
    fluxes = zeros(length(simPoints), length(model.rxns));
    
    %Load dummy data (avg nbs data with noise) for one infant
    NBSfluxdata=integrateNBSMetaboliteConcentrations(T(1,:));
    
    counter=1;
    for i = simSettings.start:simSettings.end
        
        % create infant model
        model=modelOri;
        age=i;
        fprintf(strcat('Starting infant model at age_ ',string(age), ' with weight_',string(weight) ))
        model=create_newborn(model, weight, age, NBSfluxdata);
        
        % run simulation
        lpSolution = optimizeWBModel(model)        
        
        % in case of infeasible solutions STOP the evaluation
        if lpSolution.stat == 0
            return
        end
        
        % growth equals solution of biomass reaction
        % Update weight for male and female models
        if strcmp(model.sex, 'male')
            growth=lpSolution.f;
            fluxes(counter,:) = lpSolution.v;
            lpSolution.weight=weight;
            weight = weight* growth; 
        else
            growth=round(lpSolution.f,4);
            fluxes(counter,:) = lpSolution.v;
            lpSolution.weight=weight;
            weight = round(weight* growth,4); 
        end
        
        weights(counter) = weight*1000;
        counter=counter+1;
        
        %Save intermediate results and models
        
        simulationResults = storeResults(simPoints, weights, fluxes, model);
        %resultsPATH=strcat(currentfolder,'/results/3_2_Results/'); %USER INPUT
        %save(strcat(resultsPATH, model.sex,'_growth_models/results/growth_checkpoint.mat'), 'simulationResults'); 
        %Store LP solution at every time step is necessary
        %save(strcat(resultsPATH, model.sex,'_growth_models/models/',model.sex,'_growth_day',int2str(age),'.mat'), 'lpSolution');

    end
    
    simulationResults = storeResults(simPoints, weights, fluxes, model);
end



function simPoints = getSimPoints(simSettings)
    startPoint = simSettings.start;
    endPoint = simSettings.end;
    simPoints = startPoint:endPoint;
end


function simulationResults = storeResults(simPoints, weights, fluxes, model)
    simulationResults.timePoints = simPoints;   
    simulationResults.weight = weights;
    simulationResults.v = fluxes;
    simulationResults.rxns=model.rxns;
end


