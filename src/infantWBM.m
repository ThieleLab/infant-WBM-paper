function model_new=infantWBM(sex,age)
%function model=infantWBM(sex,age)
% This function creates a infantWBM based on a given gender and age.
% It loads the infant's weight corresponding to the age and calls the
% create_newborn function to create an infant model. 

%INPUT
% sex       gender of the model either 'male' or 'female'
% age       age of the infant in days, e.g. '3' 

%OUTPUT
% model_new     personalized infant model

% Elaine Zaunseder 2023

% check gender of the model and load corresponding data 
    if strcmp(sex,'male')
        load('models/babyBoy.mat')
        initialWeight=3.3;
        weights=readtable('/results/3_2_Results/male_growth_models/reference_data/baby_boy_weight.csv');
        T= readtable('dummy_normal_male.xlsx');
    elseif strcmp(sex, 'female')
        load('models/babyGirl.mat')
        initialWeight=3.2;
        weights=readtable('/results/3_2_Results/female_growth_models/reference_data/baby_girl_weight.csv');
        T= readtable('dummy_normal_female.xlsx');
    else
        fprintf('Please enter a valid gender, these are "male" or "female"!')
        return
    end
    
% Load input data for models (metabolite concentrations)
    NBSfluxdata=integrateNBSMetaboliteConcentrations(T(1,:));
    
% Load weight corresponding to the age
    if age<=1
        weight=initialWeight;
    elseif age <=180
        weight=round(weights{age-1,1}/1000,4);
    else
        fprintf('Please provide an age in the range 1-180 days!')
        return
    end
    
% Create infant-WBM
    model_new=create_newborn(model, weight, age, NBSfluxdata);

end