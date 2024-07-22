function [IEMSolutions,IEMTable]= timepointIEMAnalysis(person,resultsPATH,timepoints, sex, geneMarkerList, urine, compartment)
% This function applies the performIEMAnalysis.m for several time points and IEMs
% (given in the geneMarkerList) for an infant-WBM or adult WBM and saves 
% the results such as the IEMTable in a given results path. 

% INPUT
% person            Determines, whether the model is an infant model 'infant'
%                   or an adult model 'adult'
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

% Author: Elaine Zaunseder October 2023

% IEM Analysis for several time points 
    for i=1:length(timepoints)
        point=1;
        if strcmp(person,'infant')
            
            %load model
            %load(strcat(sex,'_day_',string(timepoints(i)),'.mat'))
            % create model
            model=infantWBM(sex, i);

            % Set constraints
            Whole_ID = find(contains(model.rxns,'Whole')==1);
            model = changeObjective(model,model.rxns(Whole_ID));
            
            %make feasibility adjustments
            model.lb(ismember(model.rxns,'Diet_EX_phe_L[d]'))=-10;
            if timepoints(i) == 1
                model.lb(ismember(model.rxns,'Diet_EX_tyr_L[d]'))=-1.5;
            end
         
           if strcmp(sex,'male')
                %make feasibility adjustments for male models
                model.lb(ismember(model.rxns,'Diet_EX_xyl_D[d]'))=-1;
                
                %optimize WT model and set constraints
                solution = optimizeCbModel(model,'max');
                model.lb(Whole_ID)=1;%solution.f;
                model.ub(Whole_ID)=solution.f;
                
           elseif strcmp(sex, 'female')
                solution = optimizeCbModel(model,'max');
                
                %optimize WT model and set constraints
                model.lb(Whole_ID)=0.5;%solution.f;
                model.ub(Whole_ID)=solution.f;
           end          
            
        %%%%% TEST CODE NOT ACTIVE %%%%   
        elseif strcmp(person,'nadult') %Test code for adults, not active at the moment
            if strcmp(sex,'female')
               % load('Harvetta_104b.mat');  EZ adult
                load Harvetta_1_03c.mat
                model=female;
            elseif strcmp(sex,'male')
               % load('Harvey_104b.mat'); EZ adult
                load Harvey_1_04c.mat
                model=male;
            end
            Whole_ID = find(contains(model.rxns,'Whole')==1); %EZ adult
           % model = changeObjective(model,model.rxns(Whole_ID)); EZ adult
          %  model.lb(Whole_ID)=1; % EZ adult
           % model.ub(Whole_ID)=1; % EZ adult
           
           %%% Start Ines code for adults (IEM_biomarker_day_simulation.mlx)
           %% set unified reaction constraints -- they are duplicated again in individual scripts

            R = {'_ARGSL';'_GACMTRc';'_FUM';'_FUMm';'_HMR_7698';'_UAG4E';'_UDPG4E';'_GALT'; '_G6PDH2c';'_G6PDH2r';'_G6PDH2rer';...
                '_GLUTCOADHm';'_r0541'; '_ACOAD8m';'_RE2410C';'_RE2410N'};
            RxnsAll2 = '';
            for i = 1: length(R)
                RxnsAll = model.rxns(~cellfun(@isempty,strfind(model.rxns,R{i})));
                RxnsAll2 =[RxnsAll2;RxnsAll];
            end

            %excluded reactions
            R2 = {'_FUMt';'_FUMAC';'_FUMS';'BBB'};
            RxnsAll4 = '';
            for i = 1: length(R2)
                RxnsAll3 = model.rxns(~cellfun(@isempty,strfind(model.rxns,R2{i})));
                RxnsAll4 =[RxnsAll4;RxnsAll3];
            end
            RxnsAll4 = unique(RxnsAll4);
            IEMRxns = setdiff(RxnsAll2,RxnsAll4);
            RxnMic = model.rxns(~cellfun(@isempty,strfind(model.rxns,'Micro_')));
            if ~isempty(RxnMic)
                RxnMic
            end
            IEMRxns = setdiff(IEMRxns,RxnMic);
            % set ARGSL to be irreversible
            model.lb(ismember(model.rxns,IEMRxns)) = 0;

            R2 = {'_r0784';'_r0463'};
            RxnsAll2 = '';
            for i = 1: length(R2)
                RxnsAll = model.rxns(~cellfun(@isempty,strfind(model.rxns,R2{i})));
                RxnsAll2 =[RxnsAll2;RxnsAll];
            end
            X = unique(RxnsAll2);
            RxnMic = model.rxns(~cellfun(@isempty,strfind(model.rxns,'Micro_')));
            if ~isempty(RxnMic)
                RxnMic
            end
            X = setdiff(X,RxnMic);
            model.lb(ismember(model.rxns,X)) = 0;
            model.ub(ismember(model.rxns,X)) = 0;

            Rnew = {'BileDuct_EX_12dhchol[bd]_[luSI]';'BileDuct_EX_3dhcdchol[bd]_[luSI]';'BileDuct_EX_3dhchol[bd]_[luSI]';'BileDuct_EX_3dhdchol[bd]_[luSI]';'BileDuct_EX_3dhlchol[bd]_[luSI]';'BileDuct_EX_7dhcdchol[bd]_[luSI]';'BileDuct_EX_7dhchol[bd]_[luSI]';'BileDuct_EX_cdca24g[bd]_[luSI]';'BileDuct_EX_cdca3g[bd]_[luSI]';'BileDuct_EX_cholate[bd]_[luSI]';'BileDuct_EX_dca24g[bd]_[luSI]';'BileDuct_EX_dca3g[bd]_[luSI]';'BileDuct_EX_dchac[bd]_[luSI]';'BileDuct_EX_dgchol[bd]_[luSI]';'BileDuct_EX_gchola[bd]_[luSI]';'BileDuct_EX_hca24g[bd]_[luSI]';'BileDuct_EX_hca6g[bd]_[luSI]';'BileDuct_EX_hdca24g[bd]_[luSI]';'BileDuct_EX_hdca6g[bd]_[luSI]';'BileDuct_EX_hyochol[bd]_[luSI]';'BileDuct_EX_icdchol[bd]_[luSI]';'BileDuct_EX_isochol[bd]_[luSI]';'BileDuct_EX_lca24g[bd]_[luSI]';'BileDuct_EX_tchola[bd]_[luSI]';'BileDuct_EX_tdchola[bd]_[luSI]';'BileDuct_EX_tdechola[bd]_[luSI]';'BileDuct_EX_thyochol[bd]_[luSI]';'BileDuct_EX_uchol[bd]_[luSI]'};
            model.ub(ismember(model.rxns,Rnew)) = 100;


            model = changeObjective(model,model.rxns(Whole_ID));
            model.lb(Whole_ID)=1;
            model.ub(Whole_ID)=1;
            
            %load model
            load(strcat(sex,'_day_',string(timepoints(i)),'.mat'))
            Whole_ID = find(contains(model.rxns,'Whole')==1);
            model = changeObjective(model,model.rxns(Whole_ID));
            
            %make feasibility adjustments
            model.lb(ismember(model.rxns,'Diet_EX_phe_L[d]'))=-10;
            if timepoints(i) == 1
                model.lb(ismember(model.rxns,'Diet_EX_tyr_L[d]'))=-1.5;
            end
         
           if strcmp(sex,'male')
                %make feasibility adjustments for male models
                model.lb(ismember(model.rxns,'Diet_EX_xyl_D[d]'))=-1;
                
                %optimize WT model and set constraints
                solution = optimizeCbModel(model,'max');
                model.lb(Whole_ID)=1;%solution.f;
                model.ub(Whole_ID)=solution.f;
                
           elseif strcmp(sex, 'female')
                solution = optimizeCbModel(model,'max');
                
                %optimize WT model and set constraints
                model.lb(Whole_ID)=0.5;%solution.f;
                model.ub(Whole_ID)=solution.f;
           end          
            
        %%% CODE FOR ADULT MODELS %%%%   
        elseif strcmp(person,'adult')
            if strcmp(sex,'female')
               % load('Harvetta_104b.mat');  EZ adult
                load Harvetta_1_03c.mat
                model=female;
            elseif strcmp(sex,'male')
               % load('Harvey_104b.mat'); EZ adult
                load Harvey_1_04c.mat
                model=male;
            end
             Whole_ID = find(contains(model.rxns,'Whole')==1); %EZ adult
           % model = changeObjective(model,model.rxns(Whole_ID)); EZ adult
          %  model.lb(Whole_ID)=1; % EZ adult
           % model.ub(Whole_ID)=1; % EZ adult
           
           %%% Start Ines code for adults (IEM_biomarker_day_simulation.mlx)
           %% set unified reaction constraints -- they are duplicated again in individual scripts

            R = {'_ARGSL';'_GACMTRc';'_FUM';'_FUMm';'_HMR_7698';'_UAG4E';'_UDPG4E';'_GALT'; '_G6PDH2c';'_G6PDH2r';'_G6PDH2rer';...
                '_GLUTCOADHm';'_r0541'; '_ACOAD8m';'_RE2410C';'_RE2410N'};
            RxnsAll2 = '';
            for i = 1: length(R)
                RxnsAll = model.rxns(~cellfun(@isempty,strfind(model.rxns,R{i})));
                RxnsAll2 =[RxnsAll2;RxnsAll];
            end

            %excluded reactions
            R2 = {'_FUMt';'_FUMAC';'_FUMS';'BBB'};
            RxnsAll4 = '';
            for i = 1: length(R2)
                RxnsAll3 = model.rxns(~cellfun(@isempty,strfind(model.rxns,R2{i})));
                RxnsAll4 =[RxnsAll4;RxnsAll3];
            end
            RxnsAll4 = unique(RxnsAll4);
            IEMRxns = setdiff(RxnsAll2,RxnsAll4);
            RxnMic = model.rxns(~cellfun(@isempty,strfind(model.rxns,'Micro_')));
            if ~isempty(RxnMic)
                RxnMic
            end
            IEMRxns = setdiff(IEMRxns,RxnMic);
            % set ARGSL to be irreversible
            model.lb(ismember(model.rxns,IEMRxns)) = 0;

            R2 = {'_r0784';'_r0463'};
            RxnsAll2 = '';
            for i = 1: length(R2)
                RxnsAll = model.rxns(~cellfun(@isempty,strfind(model.rxns,R2{i})));
                RxnsAll2 =[RxnsAll2;RxnsAll];
            end
            X = unique(RxnsAll2);
            RxnMic = model.rxns(~cellfun(@isempty,strfind(model.rxns,'Micro_')));
            if ~isempty(RxnMic)
                RxnMic
            end
            X = setdiff(X,RxnMic);
            model.lb(ismember(model.rxns,X)) = 0;
            model.ub(ismember(model.rxns,X)) = 0;

            Rnew = {'BileDuct_EX_12dhchol[bd]_[luSI]';'BileDuct_EX_3dhcdchol[bd]_[luSI]';'BileDuct_EX_3dhchol[bd]_[luSI]';'BileDuct_EX_3dhdchol[bd]_[luSI]';'BileDuct_EX_3dhlchol[bd]_[luSI]';'BileDuct_EX_7dhcdchol[bd]_[luSI]';'BileDuct_EX_7dhchol[bd]_[luSI]';'BileDuct_EX_cdca24g[bd]_[luSI]';'BileDuct_EX_cdca3g[bd]_[luSI]';'BileDuct_EX_cholate[bd]_[luSI]';'BileDuct_EX_dca24g[bd]_[luSI]';'BileDuct_EX_dca3g[bd]_[luSI]';'BileDuct_EX_dchac[bd]_[luSI]';'BileDuct_EX_dgchol[bd]_[luSI]';'BileDuct_EX_gchola[bd]_[luSI]';'BileDuct_EX_hca24g[bd]_[luSI]';'BileDuct_EX_hca6g[bd]_[luSI]';'BileDuct_EX_hdca24g[bd]_[luSI]';'BileDuct_EX_hdca6g[bd]_[luSI]';'BileDuct_EX_hyochol[bd]_[luSI]';'BileDuct_EX_icdchol[bd]_[luSI]';'BileDuct_EX_isochol[bd]_[luSI]';'BileDuct_EX_lca24g[bd]_[luSI]';'BileDuct_EX_tchola[bd]_[luSI]';'BileDuct_EX_tdchola[bd]_[luSI]';'BileDuct_EX_tdechola[bd]_[luSI]';'BileDuct_EX_thyochol[bd]_[luSI]';'BileDuct_EX_uchol[bd]_[luSI]'};
            model.ub(ismember(model.rxns,Rnew)) = 100;


            model = changeObjective(model,model.rxns(Whole_ID));
            minRxnsFluxHealthy=0.75;
            clear IEMSolutions IEMTable missingMet
            causal = 1;
            urine = 0;
            model.lb(contains(model.rxns,'HMR_0761')) = 0;
            model.ub(contains(model.rxns,'HMR_0761')) = 0;
            %% end Ines code

        end

        % run IEM analysis
        minRxnsFluxHealthy=0.75;
        clear IEMSolutions IEMTable missingMet
        causal = 1;
        
        [IEMSolutions,IEMTable] = performIEMAnalysis(model,geneMarkerList,compartment,urine,minRxnsFluxHealthy,causal,'ibm_cplex',timepoints(point));

        % save results
       % save(strcat(resultsPATH,'IEMTable',string(timepoints(i)),'_',person,sex),'IEMTable')
       % save(strcat(resultsPATH,'IEMSolutions',string(timepoints(i)),'_',person,sex),'IEMSolutions')
    end

end