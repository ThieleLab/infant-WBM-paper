function Metabolome=FASTcomputeIndividualInfants (modelOri, age, T, start, end_num, resultsPATH)
% This function computes the LP of several individuals for specific growth
% rates given in the growthArray. This is faster because the whole-body
% biomass reaction is set to a fixed value in this case. 
%
% function Metabolome=FASTcomputeIndividualInfants (modelOri, age, T, start, number_of_models, resultsPATH)
%
% INPUT
% modelOri      model structure (boy or girl infant-WBM) 
% age           age in days
% T             table with newborn screening data from several individuals
% start         number of model to start with
% end_num       number of model to stop evaluation
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
    WT_fluxes.v=zeros(length(modelOri.rxns), end_num-start+1);
    WT_fluxes.weight=T.GEWICHT(start:end_num);
    WT_fluxes.number=T.Var1(start:end_num);
    counter=1;
    
    %establich growth rate array for growth rates to be tested. 
    growthArray={1.0092; 1.0091; 1.009; 1.0089};
    
    %start evaluating models
    for i= start:end_num
        model=modelOri;

        % Individual blood concentration data and weight from NBS
        NBSfluxdata=integrateNBSMetaboliteConcentrations(T(i,:));
        weight=T.GEWICHT(i)/1000;
        
        % set NBSflux input fix to reduce computation time
        NBSfluxdata(:,2)=NBSfluxdata(:,3);

        % Create newborn
        model=create_newborn(model, weight, age, NBSfluxdata);
        
        % Iterate through growth array and check whether model is feasible
        for k=1:length(growthArray)
            model.lb(contains(model.rxns,'Whole'))=growthArray{k};
            model.ub(contains(model.rxns,'Whole'))=growthArray{k};
            % Compute Flux 
            [solution_model] = optimizeWBModel(model);
            if solution_model.origStat==3
                continue
            else
                break
            end
        end
  
        % Check status of model solution and save results
        if solution_model.origStat==3
            WT_fluxes.v(:,counter)=zeros(length(model.rxns),1)+0.5;
        elseif solution_model.origStat==11
            WT_fluxes.v(:,counter)=zeros(length(model.rxns),1);
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
            save(strcat(resultsPATH,'Individualscheckpoint_',string(start)), 'Metabolome', '-v7.3')
            fprintf('\n Finished model %s\n', string(i));
        end
    end

    % Save final results 
    Metabolome=[WT_fluxes];
    %save(strcat(resultsPATH, 'Individuals_',string(start), '_', string(end_num)), 'Metabolome', '-v7.3')

end