function model = NewbornsetDietConstraints(model, Diet, factor, milk, age)
% This function sets diet constraints onto the bounds of the diet uptake
% reactions of a whole-body metabolic model for newborns and infants.
% Units are given in mmol/day/person.
%
% function model = setDietConstraints(model, Diet, factor, milk, age)
%
% INPUT
% model         model structure 
% Diet          diet as created by getMetaboliteFlux
% factor        value between 0 and 1; default is 1, i.e, 100% of the provided diet
% milk          milk intake
% age           age in days of infant
%
% OUTPUT
% model         updated model struct
% 
% Ines Thiele 2016-2019
% Elaine Zaunseder 2022 (update for infants)

if ~exist('factor','var')
    factor = 1; % 100% of the provided diet
end

%Adaptations for feasibility of male and female infant diet
newborn_diet=Diet;
% identified diet components and boundaries
diet_components= {'Diet_EX_chol[d]', 'Diet_EX_met_L[d]', 'Diet_EX_pe_hs[d]', 'Diet_EX_hmcr[d]', ...
    'Diet_EX_glc_D[d]', 'Diet_EX_ile_L[d]',  'Diet_EX_lys_L[d]', 'Diet_EX_val_L[d]', ...
    'Diet_EX_phe_L[d]', 'Diet_EX_thr_L[d]',  'Diet_EX_leu_L[d]', 'Diet_EX_thmmp[d]', 'Diet_EX_gudac[d]' };

% Diet bounds of these metabolite differ for female and male infants
if strcmp(model.sex,'male')
    %male infants
    diet_bounds={3.5, 8.3, 30, 0.84, 8.3, 2, 3.1, 2.4,3,2.4,5, 0.84, 0.84};
else
    %female infants
    diet_bounds={3.96, 8.3, 30, 1.18, 8.3, 2, 3.1, 2.92,3,2.4,5, 0.84, 0.84};
end

%iteratively add bounds for diet components by adding new diet components or updating old bounds
for i=1:length(diet_components)
    if ismember(newborn_diet(:,1), diet_components{i})==0
        % add new diet components
        newborn_diet{length(newborn_diet)+1,1}=diet_components{i};
        newborn_diet{length(newborn_diet),2}=diet_bounds{i};
    else
    % update old diet components
    newborn_diet{strcmp(newborn_diet(:,1), diet_components{i}),2}=max(diet_bounds{i},newborn_diet{strcmp(newborn_diet(:,1), diet_components{i}),2} );
    end
end

Diet=newborn_diet;


microbiotaEnabling =1; % consider essential microbial metabolites to be added to diet

% load essential metabolite list for AGORA models
AGORAEssentialMetabolites;
AGORAessential = regexprep(AGORAessential,'EX_','Diet_EX_');
AGORAessential = regexprep(AGORAessential,'\[u\]','\[d\]');

% set all  uptakes to 0
tmp = strmatch('Diet_EX_',model.rxns);
modelO = model;
model.lb(tmp(1:end))=0;
model.ub(tmp(1:end))=0;

% ensure uptake of metabolites required for microbiota (may not all be
% needed) -- I should really check which ones of those are need
if microbiotaEnabling == 1
    MissingUptakes = setdiff(AGORAessential,Diet(:,1));
    % open uptake for those reactions
    % set constraints to a default
    model.lb(ismember(model.rxns,MissingUptakes))=-0.1;
end

% these compounds are in the diet and needed for the proper function of the
% model but are not reported in our diet database
if 1
    MissingDietCompounds={'Diet_EX_asn_L[d]'; 'Diet_EX_gln_L[d]';'Diet_EX_chol[d]';'Diet_EX_crn[d]';'Diet_EX_elaid[d]';...
        'Diet_EX_hdcea[d]';'Diet_EX_dlnlcg[d]';'Diet_EX_adrn[d]';'Diet_EX_hco3[d]';...
        %June 2nd 2017 debug for increased brain atp
        % adding these metabolites significantly incresed the DM_atp of brain
        % -- it would be very interesting to perform a sensitivity analysis and
        % see which metabolites have which effects --> could correlate with
        % brain activity and cognition/neurodegeneration
        % many of these are also microbially produced and should have a
        % positive effect.
        % DM_atp increased from 3600 to 5015
        'Diet_EX_sprm[d]'; 'Diet_EX_carn[d]';'Diet_EX_7thf[d]';...
        'Diet_EX_Lcystin[d]';%??
        'Diet_EX_hista[d]';'Diet_EX_orn[d]';...
        'Diet_EX_ptrc[d]';'Diet_EX_creat[d]';
        };
    
    model.lb(ismember(model.rxns,MissingDietCompounds))=-50;
    MissingDietCompounds={
        % 'Diet_EX_uri[d]'
        'Diet_EX_cytd[d]'
        % 'Diet_EX_gam[d]'
        %'Diet_EX_gal[d]'
        'Diet_EX_so4[d]'
        %'Diet_EX_fuc_L[d]'
        %'Colon_DM_Asn_X_Ser_Thr_ly_'
        %'Diet_EX_h[d]'
        
        };
    model.lb(ismember(model.rxns,MissingDietCompounds))=-50;

%if strcmp(model.sex,'male')
%    model.lb(ismember(model.rxns,'Diet_EX_chol[d]'))=-41.251; %based on a daily intake of 396 mg in Av Am Diet per day (Sahoo 2013 paper)
%    'CHOL'
%end

end

% micronutrient - defined to have mole/day/person rate below 1e-6 mol/day/person
% lower bounds will be relaxed by factor 10
micronutrients ={%'Diet_EX_adpcbl[d]'
    % I changed the exchange ID since it was generated based on the
    % metabolite (ID adocbl). AH 16/12/01
    'Diet_EX_adocbl[d]'
    'Diet_EX_vitd2[d]'
    'Diet_EX_vitd3[d]'
    'Diet_EX_psyl[d]'
    'Diet_EX_gum[d]'
    'Diet_EX_bglc[d]'
    'Diet_EX_phyQ[d]'
    'Diet_EX_fol[d]'
    'Diet_EX_5mthf[d]'
    'Diet_EX_q10[d]'
    'Diet_EX_retinol_9_cis[d]'
    'Diet_EX_pydxn[d]'
    'Diet_EX_pydam[d]'
    'Diet_EX_pydx[d]'
    'Diet_EX_pheme[d]'
    'Diet_EX_ribflv[d]'
    'Diet_EX_thm[d]'
    % added 24.08.2016
    'Diet_EX_avite1[d]'
    'Diet_EX_pnto_R[d]'
    };
%no uptake enforced but max uptake rate
ions={    'Diet_EX_na1[d]'	%'0.056546816'
    'Diet_EX_cl[d]'	%'0.056412715'
    'Diet_EX_k[d]'	%'0.12020983'
   % 'Diet_EX_pi[d]'	%'0.007143207'
    'Diet_EX_zn2[d]'
    'Diet_EX_cu2[d]'
    };
so4={
    %02.06.2017
    'Diet_EX_so4[d]' % see text above about addition of met for brain demand
    };

for i = 1:  size(Diet,1)
    R = Diet(i,1);
    % exception for micronutrients to avoid numberical issues
    if ~isempty(find(ismember(micronutrients,R)))&& cell2mat(Diet(i,2))<=10
        model.lb(find(ismember(model.rxns,R))) = -10*factor;%-1.2*str2num(Diet{i,2})*factor;
    elseif ~isempty(find(ismember(micronutrients,R))) && str2num>0.1
        model.lb(find(ismember(model.rxns,R))) = -1.2*cell2mat(Diet(i,2))*100*factor;
    elseif ~isempty(find(ismember(ions,R)))
        %  model.lb(find(ismember(model.rxns,R))) = -1.2*str2num(Diet{i,2});
        model.lb(find(ismember(model.rxns,R))) = -1.2*cell2mat(Diet(i,2))*100*factor;
    elseif ~isempty(find(ismember(so4,R)))
        model.lb(find(ismember(model.rxns,R))) = -1000*factor;
    else
        model.lb(find(ismember(model.rxns,R))) = -1.2*cell2mat(Diet(i,2))*factor;
    end
end
% do not enforce diet uptake --> set if statement to 0
if 1
    for i = 1:size(Diet,1) % fine until 70
        R = Diet{i,1};
        if ~isempty(find(ismember(micronutrients,R)))
            model.ub(find(ismember(model.rxns,R))) =  -0.8*cell2mat(Diet(i,2))*factor;
        elseif ~isempty(find(ismember(ions,R)))
            model.ub(find(ismember(model.rxns,R))) = -0.8*cell2mat(Diet(i,2))*factor;
        else
            model.ub(find(ismember(model.rxns,R))) = -0.8*cell2mat(Diet(i,2))*factor;
        end
    end
end

% Adaptations for feasibility of solution
model.lb(ismember(model.rxns,'Diet_EX_5mthf[d]'))=-0.1; 
model.lb(ismember(model.rxns,'Diet_EX_phyQ[d]'))=-0.1; 
model.lb(ismember(model.rxns,'Diet_EX_thf[d]'))=-0.1;
model.lb(ismember(model.rxns,'Diet_EX_pi[d]'))=-214;

% Load data from table for Lysine and Lcystine diet intake
T=readtable('LysLcyDiet');
model.lb(strcmp(model.rxns,'Diet_EX_lys_L[d]'))=T{age,'Lysine'};
%Lcystine
if strcmp(model.sex,'female')
    model.lb(ismember(model.rxns,'Diet_EX_Lcystin[d]'))=T{age,'Lcystin'};
end

%Initial milk diet did not contain lactose, hence, we add 7% of lactose 
milk_intake=milk(age,2);
lactose=((milk_intake*0.07)/342.3)*1000;
model.lb(ismember(model.rxns,'Diet_EX_lcts[d]'))=-lactose;


model.SetupInfo.DietComposition = Diet;

