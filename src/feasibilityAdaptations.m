function model=feasibilityAdaptations(model)
% For the infant-WBM models some reactions were identified which limited
% the feasibility of the models.
%These constraints are adapted in this script to enable feasibility for all
%models. 

%Ines Thiele & Elaine Zaunseder 2023

model.lb(contains(model.rxns,'EX_nh4[u]'))=0;

model.lb(contains(model.rxns,'EX_co2[a]')) = 500;
model.ub(contains(model.rxns,'EX_o2[a]')) = -4000; % chosen based on FBA solution feasible

model.lb(contains(model.rxns,'Muscle_EX_ala_L(e)_[bc]'))=0;

% somehow the physiologicalDefaultP script does not update the ub
model.ub(contains(model.rxns,'Muscle_EX_o2(e)_[bc]')) = 0;

model.lb(contains(model.rxns,'RBC_EX_dag_hs(e)_[bc]')) = -1;
model.lb(contains(model.rxns,'RBC_EX_h2o2(e)_[bc]')) = 0;
model.ub(contains(model.rxns,'RBC_EX_h2o2(e)_[bc]')) = 0;
model.lb(contains(model.rxns,'RBC_EX_hco3(e)_[bc]')) = 0;
model.lb(contains(model.rxns,'RBC_EX_co(e)_[bc]')) = 0;

model.ub(contains(model.rxns,  'Brain_EX_o2(e)_[csf]')) = 0;


end