function [QP]= QPInfants(gender)


param.minNorm = 1e-6;
timepoints=[1 30 60 90 120 150 180]; % set time points for every month (e.g. every 30 days)

load("QPgrowth.mat")

% for sex-specific models and growth rates
if strcmp(gender,'male')
    load("babyBoy.mat")   
    QPgrowthrate=QPgrowth.male;
    param.minNorm = 1e-6;
else
    load('babyGirl.mat')
    QPgrowthrate=QPgrowth.female;
    param.minNorm = 1e-6;
end    

QP.rxns=model.rxns;
QP.v=zeros(180,length(model.rxns));

%run QP for several time steps
for i=timepoints
    fprintf(strcat('\n Model day_',string(i),'\n'))
    %load(strcat(gender,'_DAY_',string(i),'.mat'))
    model=infantWBM(gender,i);
    model.lb(contains(model.rxns,'Whole'))=QPgrowthrate(i);
    model.ub(contains(model.rxns,'Whole'))=QPgrowthrate(i);
    save(strcat('Female_day_',string(i),'.mat'), 'model')
   % qpSolution = optimizeWBModel_QP(model, param)
    qpSolution = optimizeWBModel(model, param)
    QP.v(i,:) = qpSolution.v;
end
  

timepoints=[1 30 60 90 120 150 180]; % set time points for every month (e.g. every 30 days)

load("QPgrowth.mat")

% for sex-specific models and growth rates
if strcmp(gender,'male')
    load("babyBoy.mat")   
    QPgrowthrate=QPgrowth.male;
    param.minNorm = 1e-6;
else
    load('babyGirl.mat')
    QPgrowthrate=QPgrowth.female;
    param.minNorm = 1e-6;
end    

QP.rxns=model.rxns;
QP.v=zeros(180,length(model.rxns));

%run QP for several time steps
for i=timepoints
    fprintf(strcat('\n Model day_',string(i),'\n'))
    %load(strcat(gender,'_DAY_',string(i),'.mat'))
    model=infantWBM(gender,i);
    model.lb(contains(model.rxns,'Whole'))=QPgrowthrate(i);
    model.ub(contains(model.rxns,'Whole'))=QPgrowthrate(i);
    save(strcat('Female_day_',string(i),'.mat'), 'model')
   % qpSolution = optimizeWBModel_QP(model, param)
    qpSolution = optimizeWBModel(model, param)
    QP.v(i,:) = qpSolution.v;
end
   
end