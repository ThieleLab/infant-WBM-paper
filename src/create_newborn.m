function [model] = create_newborn(model, weight, age, NBSfluxdata)
%function model = create_newborn(model, weight, age, NBSfluxdata)
%function to create age, weight and sex dependent newborn WBM model

%INPUT
% model         a WBM model
% weight        weight in kg
% age           age in days
% NBSfluxdata   Blood concentration data from a newborn

%OUTPUT
% model      constrained model representing newborn

% Elaine Zaunseder Sep 2022 - May 2023

%check if weight is given in kg
if weight > 100
    fprintf('ERROR: Weight input has to be given in kg!')
    return
end

% 1. INDIVIDUAL PARAMETERS
IndividualParameters=newbornPhysiolDefaultParameters(model.sex, weight, age); 
IndividualParameters.OrgansWeights(:,2)=num2cell(cell2mat(IndividualParameters.OrgansWeights(:,3))*weight*1000); % organ growth depending on weight

%set physiological constraints
model = physiologicalConstraintsHMDBbased(model,IndividualParameters); %takes ~25 seconds
model = physiologicalConstraintsHMDBbased(model,IndividualParameters, [], 'direct', NBSfluxdata, 'bc'); %takes ~3.5 seconds


% 2. DAILY DIET takes ~1.5 seconds

milk=load('milkModel_nbs.txt');  %milk input from STIG MET Paper https://www.nature.com/articles/s41540-017-0004-5
milk_intake=milk(age,2);
intake_percentage=milk_intake/233;
newborn_diet=getMetaboliteFlux({'Food_EX_Milk_Human_Mature_Fluid[d]'  milk_intake/100}); % VMH breast milk diet per 100 ml
model=NewbornsetDietConstraints(model, newborn_diet,1, milk, age); % set diet constraints in model

% Assuming newborn is excreting 10% of water from adult at age day 2, where milk intake is 233 ml
model = changeRxnBounds(model,'EX_h2o[a]',4718*intake_percentage*0.9,'l'); %air
model = changeRxnBounds(model,'EX_h2o[a]',4718*intake_percentage*1.1,'u');
model = changeRxnBounds(model,'EX_h2o[sw]',3608*intake_percentage*0.9,'l'); %sweat
model = changeRxnBounds(model,'EX_h2o[sw]',3608*intake_percentage*1.1,'u');
model = changeRxnBounds(model,'EX_h2o[u]',7771*intake_percentage*0.9,'l'); %urine
model = changeRxnBounds(model,'EX_h2o[u]',7771*intake_percentage*1.1,'u');
model = changeRxnBounds(model,'Excretion_EX_h2o[fe]',555*intake_percentage*0.9,'l'); % feces
model = changeRxnBounds(model,'Excretion_EX_h2o[fe]',555*intake_percentage*1.1,'u');


% 3. DAILY ATP DEMAND
% Muscles
activityModel = makeActivityModel('inputs/activityModel.txt', 14.4); % activity model from STIG MET Paper https://www.nature.com/articles/s41540-017-0004-5
atp_activity=((activityModel(age)*weight)/28.1250)*1000;
model.lb(ismember(model.rxns,'Muscle_DM_atp_c_'))=atp_activity;

% Brain, Heart, Adipocytes(Thermoregulation)
heart_weight=IndividualParameters.OrgansWeights{contains(IndividualParameters.OrgansWeights(:,1),'Heart'),2};
heart_percentage=heart_weight/331; % adult male heart 331 gram heart with bound 6000 in adult model
model = changeRxnBounds(model,'Heart_DM_atp_c_',6000* heart_percentage,'l'); 
brainATP=readtable('brainATPdemand');
model = changeRxnBounds(model,'Brain_DM_atp_c_',brainATP{age,2},'l'); % brain atp starting from 3000mmol  at age 0 % https://pubmed.ncbi.nlm.nih.gov/25157149/
model = changeRxnBounds(model,'Adipocytes_DM_atp_c_',brainATP{age,2}*0.33,'l');


% 4. BIOMASS REACTION
% adjust biomass coefficients depending on organ weights
model=adjustWholeBodyRxnCoeff(model, IndividualParameters.OrgansWeights(:,1), cell2mat(IndividualParameters.OrgansWeights(:,3)));
model.lb(contains(model.rxns,'Whole'))=1;
model.ub(contains(model.rxns,'Whole'))=2;
model.osenseStr='max';
model = changeObjective(model, 'Whole_body_objective_rxn');


% 5. MODEL DESCRIPTION
% Add description to model to save individual parameters for personalization for each model
model.status='adapted';
model.modelName='Newborn model adapted with age, sex, weight and diet';
model.diet={'milkmodel(gr)',milk(age,2)};
model.IndividualParameters=IndividualParameters;


% 6. FEASIBILITY ADAPTIONS
model=feasibilityAdaptations(model);

end