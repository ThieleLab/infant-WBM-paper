function [IEMSolutions,IEMTable]= IEMDietaryTreatment(strategy,resultsPATH,timepoints, sex, geneMarkerList, urine, compartment)
% This function applies the performIEMAnalysis.m for several time points and IEMs
% (given in the geneMarkerList) for an infant-WBM or adult WBM and saves
% the results such as the IEMTable in a given results path.

% INPUT
% strategy          gives the dietary treatment strategy
%                   .target is the diet flux that is suppose to be reduced
%                   .percent is the amount of reduction (eg. -0.25) or
%                   increase (e.g. 0.25) for a 25% decrease or increase
% resultsPATH       Path to results folder
% timepoints        One or several timepoints in days
% sex               Gender of the model 'male' or 'female'
% geneMarkerList    list of genes for knock out and corresponding
%                   biomarkers e.g.,
%                   geneMarkerListMarkerList = {
%                             '5053.1' 'trp_L;actyr;phe_L;tyr_L'
%                             '249.1' '3pg;cholp;glyc3p;ethamp'
%                             };
% urine             0 (no) or 1 (yes) if biomarker should be analysed in urine
% compartment       [bc] if biomarker in blood compartment should be
%                   analysed

% OUTPUT
% IEMSolutions      Structure containing the predictions for each gene.
%                   Metabolites that are not occurring in a biofluid, will have a 'NA' in the
%                   corresponding fields
% IEMTable          Cell array containing the predictions for each gene (same
%                   content as in IEMSolutions

% Author: Elaine Zaunseder February 2024

% IEM Analysis for several time points
if strcmp(sex, 'female')
    for j= 1: length(geneMarkerList(:,1))

       % Load model and adpt for feasibility
        load(strcat(sex,'_day_',string(timepoints),'.mat'));
        Whole_ID = find(contains(model.rxns,'Whole')==1);
        model = changeObjective(model,model.rxns(Whole_ID));
        model.lb(ismember(model.rxns,'Diet_EX_phe_L[d]'))=-10;
        model.lb(ismember(model.rxns,'Diet_EX_tyr_L[d]'))=-1.5;
        
        %Set dietary treatment strategy
        model.lb(ismember(model.rxns,strategy.target))=model.lb(ismember(model.rxns,strategy.target))*(1+strategy.percent); 
        
        % set constraints on whole body biomass reaction
        solution = optimizeCbModel(model,'max');
        model.lb(Whole_ID)=0.5;%solution.f;
        model.ub(Whole_ID)=solution.f;

        % run IEM analysis
        minRxnsFluxHealthy=0.75;
        clear IEMSolutions IEMTable missingMet
        causal = 1;
        [IEMSolutions,IEMTable] = performIEMAnalysis(model,geneMarkerList(j,:),compartment,urine,minRxnsFluxHealthy,causal,'ibm_cplex',timepoints);

        % save results
        save(strcat(resultsPATH,'Dietary_IEMTable_',sex),'IEMTable')
        save(strcat(resultsPATH,'Dietary_IEMSolutions_',sex),'IEMSolutions')
       
    end
    
%%%%%%%%%%% MALE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(sex, 'male')

     for j= 1: length(geneMarkerList(:,1))
         
        % Load correct infant model
        load(strcat(sex,'_day_',string(timepoints),'.mat'));
        
        % Adapt constraints necessary for feasibility of the models in the
        % IMD case
        Whole_ID = find(contains(model.rxns,'Whole')==1);
        model = changeObjective(model,model.rxns(Whole_ID));
        model.lb(ismember(model.rxns,'Diet_EX_phe_L[d]'))=-10;
        model.lb(ismember(model.rxns,'Diet_EX_tyr_L[d]'))=-1.5;
        model.lb(ismember(model.rxns,'Diet_EX_xyl_D[d]'))=-1;
        
        % Set dietary treatment strategy
        model.lb(ismember(model.rxns,strategy.target))=model.lb(ismember(model.rxns,strategy.target))*(1+strategy.percent); 
        
        % Set constraints on whole body biomass reaction
        solution = optimizeCbModel(model,'max');
        model.lb(Whole_ID)=1;%solution.f;
        model.ub(Whole_ID)=solution.f;

        % Run IEM analysis
        minRxnsFluxHealthy=0.75;
        clear IEMSolutions IEMTable missingMet
        causal = 1;
        [IEMSolutions,IEMTable] = performIEMAnalysis(model,geneMarkerList(j,:),compartment,urine,minRxnsFluxHealthy,causal,'ibm_cplex',timepoints);

        
        % save results
        save(strcat(resultsPATH,'DietaryIEMTable__',sex),'IEMTable')
        save(strcat(resultsPATH,'DietaryIEMSolutions_',sex),'IEMSolutions')
       
    end

end

