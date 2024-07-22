function plotQPInfants(QP, fluxType)

timepoints=[1 30 60 90 120 150 180]; % set time points for every month (e.g. every 30 days)

if strcmp(fluxType, 'water')
    % Get water excretion fluxes and convert to ml
    water_intake=abs(QP.v(timepoints,(ismember(QP.rxns,'Diet_EX_h2o[d]'))))/1000*18.02;
    water_out_a=(QP.v(timepoints,(ismember(QP.rxns,'EX_h2o[a]'))))/1000*18.02;
    water_out_sw=(QP.v(timepoints,(ismember(QP.rxns,'EX_h2o[sw]'))))/1000*18.02;
    water_out_u=(QP.v(timepoints,(ismember(QP.rxns,'EX_h2o[u]'))))/1000*18.02;
    water_out_fe=(QP.v(timepoints,(ismember(QP.rxns,'Excretion_EX_h2o[fe]'))))/1000*18.02;

    % Plot results
    plot(timepoints,water_intake,'-o', 'LineWidth', 3.0,'color', '#0072BD', 'DisplayName', 'DIET h2o')
    hold on
    plot(timepoints,water_out_a,'-o', 'LineWidth', 3.0,'color','#4DBEEE', 'DisplayName', 'EX h2o air')
    hold on
    plot(timepoints,water_out_fe,'-o', 'LineWidth', 3.0, 'color', '#7B3F00', 'DisplayName', 'EX h2o faeces')
    hold on
    plot(timepoints,water_out_u,'-o', 'LineWidth', 3.0,'color','#D6E865', 'DisplayName', 'EX h2o urine')
    hold on
    plot(timepoints,water_out_sw,'-o', 'LineWidth', 3.0,'color','#ffc4e1', 'DisplayName','EX h2o sweat')
    hold off
    set(gca,'box','off')
    xlabel('Age in days','FontSize',18)
    ylabel('Water intake/excretion in (ml/infant/day)', 'FontSize',16,'FontName', 'Arial')
    title('Water Balance in infant', 'FontSize',22,'FontName', 'Arial')
    legend('Location','northwest', 'FontSize', 13)
    
elseif strcmp(fluxType,'ATPSynthase')
    % Get reaction fluxes
    liver=QP.v(timepoints,(ismember(QP.rxns,'Liver_ATPS4m')));
    brain=QP.v(timepoints,(ismember(QP.rxns,'Brain_ATPS4m')));
    muscle=QP.v(timepoints,(ismember(QP.rxns,'Muscle_ATPS4m')));
    adipocytes=QP.v(timepoints,(ismember(QP.rxns,'Adipocytes_ATPS4m')));
    heart=QP.v(timepoints,(ismember(QP.rxns,'Heart_ATPS4m')));

    %plot results

    plot(timepoints,brain, '-o','LineWidth', 2.0,'DisplayName', 'Brain')
    hold on
    plot(timepoints,adipocytes, '-o', 'LineWidth', 2.0,'DisplayName', 'Adipose tissue')
    hold on
    plot(timepoints,heart, '-o','LineWidth', 2.0,'DisplayName', 'Heart')
    hold on
    plot(timepoints,muscle, '-o','LineWidth', 2.0,'DisplayName', 'Muscle')
    hold on
    plot(timepoints,liver,'-o', 'LineWidth', 2.0,'DisplayName', 'Liver')
    hold off
    set(gca,'box','off')
    xlabel('Age in days','FontSize',18,'FontName', 'Arial')
    ylabel('ATPS4m flux (mmol/infant/day)', 'FontSize',16,'FontName', 'Arial')
    title('ATP synthase in different organs', 'FontSize',15,'FontName', 'Arial')
    legend('Location','southeast', 'FontSize', 15)

end