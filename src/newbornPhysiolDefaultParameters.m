function IndividualParameters = newbornPhysiolDefaultParameters(sex,weight,age)

% This script creates the IndividualParameters structure which contains
% standard physiological default parameters for the reference newborn.
%INPUT
% sex           'male' or 'female'
% weight        weight in kg
% age           age in days
%
%OUTPUT
% IndividualParameters


% Ines Thiele 2016-2019
% Elaine Zaunseder 2022-2023 adapted for infants

% ORGAN PARAMETERS
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2062747/?page=5
% https://journals.lww.com/amjforensicmedicine/Abstract/2016/09000/Postmortem_Body_and_Organ_Measurements_in_Neonates.10.aspx
% https://www.rch.org.au/rchcpg/hospital_clinical_guideline_index/Neonatal_Intravenous_Fluid_Management
% Muscle Mass and Composition in Malnourished Infants and Children and Changes Seen after Recovery; P. J. REEDS, A. A. JACKSON, D. PICOU, AND N. POULTER'38'
% https://academic.oup.com/jn/article/130/9/2188/4686487

fileNameOrgan = 'OrganWeightFractAge.xlsx';
[Numbers, organWeightData] = xlsread(fileNameOrgan,'OrganWeight');
OrganNames = organWeightData(6:end,1);
if strcmp(sex,'male')
    % this is the ideal baby
    weightNewborn = Numbers(2,1);
    OrganWeightNewborn = Numbers(3:end,1);
    OrganWeightFractNewborn = Numbers(3:end,1)./Numbers(2,1);
    % now scale by input baby weight
    OrganWeightNewbornReal = OrganWeightNewborn*weight/weightNewborn;
    OrganWeightFractNewbornReal = OrganWeightFractNewborn; % fraction does not change
else
    weightNewborn = Numbers(2,2);
    OrganWeightNewborn = Numbers(3:end,2);
    OrganWeightFractNewborn = Numbers(3:end,2)./Numbers(2,2);
    OrganWeightNewbornReal = OrganWeightNewborn*weight/weightNewborn;
    OrganWeightFractNewbornReal = OrganWeightFractNewborn; % fraction does not change
end
IndividualParameters.OrgansWeights =  [OrganNames num2cell(OrganWeightNewbornReal) num2cell(OrganWeightFractNewbornReal)];

% PERSONALIZED PARAMETERS
IndividualParameters.ID = 'Infant';
IndividualParameters.bodyWeight=weight;
IndividualParameters.age=age;
IndividualParameters.sex = sex;

IndividualParameters.HeartRate =140;% beats per minute % https://academic.oup.com/jn/article/130/9/2188/4686487
IndividualParameters.StrokeVolume =max(10,(1.77+0.28)*IndividualParameters.bodyWeight); % https://pubmed.ncbi.nlm.nih.gov/3080484/
IndividualParameters.CardiacOutput = IndividualParameters.HeartRate * IndividualParameters.StrokeVolume; % in ml/min = beats/min * ml/beat

%HEMATOCRIT,
%Hematocrit_array=[62.9 61 60.5 58.7 56.6 56.1 56.6 55.5 55.2 53.7]; https://academic.oup.com/jn/article/130/9/2188/4686487
IndividualParameters.Hematocrit=0.4;% this could be continously adapted in future

% creatinine concentration in urine
IndividualParameters.MConUrCreatinineMax = 61; % mg/dL  %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4843245/
IndividualParameters.MConUrCreatinineMin = 27; % mg/dL  %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4843245/

% default maximum concetration of a metabolite in blood plasma
IndividualParameters.MConDefaultBc = 20; % uM abretary chosen

% default maximum concetration of a metabolite in csf
IndividualParameters.MConDefaultCSF = 20; % uM abretary chosen

% default maximum concetration of a metabolite in Ur
IndividualParameters.MConDefaultUrMax = 20; % umol/mmolcreatinine abretary chosen
IndividualParameters.MConDefaultUrMin = 0; % umol/mmolcreatinine

% CSF Flow rate
IndividualParameters.CSFFlowRate = (2.78-2.23 * 0 + 0.97 * log(IndividualParameters.age/365) +2.26* log(IndividualParameters.bodyWeight))/60; %ml/min % https://pubmed.ncbi.nlm.nih.gov/11818742/

% CSF to venous blood flow rate (ml/min)
IndividualParameters.CSFBloodFlowRate = 0.52;

% Urine flow rate
IndividualParameters.UrFlowRate = 2*IndividualParameters.bodyWeight*24; %ml/day, 1 – 3 ml/kg/hr --> used 2 and multiply with 24 for day
%https://www.rch.org.au/rchcpg/hospital_clinical_guideline_index/Neonatal_Intravenous_Fluid_Management/#:~:text=With%20a%20changing%20GFR%20and%20variable%20urine%20concentration%2C,offer%20additional%20information%20as%20to%20urine%20concentrating%20ability.

% GFR = Glomerular filtration rate
IndividualParameters.GlomerularFiltrationRate = 40 ;%ml/min, https://www.sciencedirect.com/science/article/pii/S1055858613000607?via%3Dihub

% Blood Flow rate per organ
fileNameOrgan = 'OrganWeightFractAge.xlsx';
[Numbers, bloodFlowData]  = xlsread(fileNameOrgan,'OrganFlowFractionScaled'); 
NumbersC = num2cell(Numbers);
bloodFlowData(3:end,2:end)=NumbersC;
exclude = {'WBC';'Lymphocytes';'Bcells';'CD4Tcells';'CD8Tcells';'Nkcells';'Monocyte';'Platelet';'RBC'};
bloodFlowData(contains(bloodFlowData(:,1),exclude),:) = [];
IndividualParameters.bloodFlowData = bloodFlowData;
IndividualParameters.bloodFlowRow = 5;
IndividualParameters.bloodFlowPercCol = [2 3];
IndividualParameters.bloodFlowOrganCol = 1;
