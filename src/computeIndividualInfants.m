function Metabolome=computeIndividualInfants (modelOri, age, T, start, number_of_models, resultsPATH)

% This function computes the LP of several individuals predicting the
% growth rate between 1 - 2. 
% This takes longer than the FASTcomputeIndividualInfants  because the whole-body
% biomass reaction is NOT set to a fixed value in this case. 
%
% function Metabolome=computeIndividualInfants (modelOri, age, T, start, number_of_models, resultsPATH)
%
% INPUT
% modelOri      model structure (boy or girl infant-WBM) 
% age           age in days
% T             table with newborn screening data from several individuals
% start         number of model to start with
% number_of_models      number of model to stop evaluation
% resultsPATH   path where results are stored
%
% OUTPUT
% Metabolome    results struct with
%      .rxns    model reactions
%      .v       solution flux vectors
%      .weight  birth weight
%      .number  ID of indidvidual data
% 
% Elaine Zaunseder 2023

    % Create structs to save results
    WT_fluxes=struct();
    WT_fluxes.rxns=modelOri.rxns;
    WT_fluxes.v=zeros(length(modelOri.rxns), number_of_models-start+1);
    WT_fluxes.weight=T.GEWICHT(start:number_of_models);
    WT_fluxes.number=T.Var1(start:number_of_models);
    infeasible=zeros(max(1,number_of_models-start),1);
    timeLimit=zeros(max(1,number_of_models-start),1);
    counter=1;
    
    %start evaluating models
    for i= start:number_of_models
        model=modelOri;

        % Individual blood concentration data and weight from NBS
        NBSfluxdata=integrateNBSMetaboliteConcentrations(T(i,:));
        NBSfluxdata(:,2)=NBSfluxdata(:,3);
        weight=T.GEWICHT(i)/1000;

        % Create newborn
        model=create_newborn(model, weight, age, NBSfluxdata);

        % Compute Flux 
        [solution_model] = optimizeWBModel(model);
        
        % Check whether problem is feasible
        if solution_model.origStat==3
            infeasible(counter)=WT_fluxes.number(counter);
            WT_fluxes.v(:,counter)=zeros(length(model.rxns),1)+0.5;
        elseif solution_model.origStat==11
            timeLimit(counter)=WT_fluxes.number(counter);
        else
            WT_fluxes.v(:,counter)=solution_model.v;
            growth=solution_model.v(contains(model.rxns, 'Whole'))
        end
        
        %increase counter
        counter=counter+1;
        
        %save checkpoint
        if ~mod(i,50)
            WT_fluxes.v(abs(WT_fluxes.v)<=0.000001)=0;
            Metabolome=[WT_fluxes];
            save(strcat(resultsPATH,'Individualscheckpoint'), 'Metabolome', '-v7.3')

            fprintf('\n Starting model %s\n', string(i+1));
        end
    end

    Metabolome=[WT_fluxes];
    save(strcat(resultsPATH, 'Individuals_', model.sex,'_',string(start), '_', string(number_of_models)), 'Metabolome', '-v7.3')

end