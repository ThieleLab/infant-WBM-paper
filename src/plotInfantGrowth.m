function plotInfantGrowth(results,gender)

%male plot
if strcmp(gender, 'male')

    STIGMet=readtable('STIGMetBoy_weight.csv');
    WHO=readtable('weightBoy.txt');

    fig = figure(); 
    fig.Position(3:4) = [700,400];
    x2 = [results.timePoints, fliplr(results.timePoints)];
    inBetween = [WHO{1:length(results.timePoints),2}', fliplr(WHO{1:length(results.timePoints),6}')];
    fill(x2, inBetween*1000,[0.96 0.96 0.96],  'DisplayName', 'WHO 0-100% quartile');
    hold on
    inBetween2 = [WHO{1:length(results.timePoints),3}', fliplr(WHO{1:length(results.timePoints),5}')];
    fill(x2, inBetween2*1000,[0.91 0.91 0.92], 'DisplayName', 'WHO 25-75% quartile');
    hold on
    plot(results.timePoints,results.weight, 'LineWidth',4, 'DisplayName', 'infantWBM')
    hold on
    plot(results.timePoints,STIGMet{1:length(results.timePoints),1}, 'LineWidth',4, 'DisplayName', 'STIGMet')
    hold off

    set(gca,'box','off')
    xlabel('Age in days','FontSize',18)
    ylabel('Body weight (g)', 'FontSize',18,'FontName', 'Arial')
    title('Male infant growth prediction', 'FontSize',22,'FontName', 'Arial')
    legend('Location','northwest', 'FontSize', 13)
    
%female plot    
else
    STIGMet=readtable('STIGMetFemale_weight.csv');
    WHO=readtable('weightFemale.txt');

    fig = figure(); 
    fig.Position(3:4) = [700,400];
    x2 = [results.timePoints, fliplr(results.timePoints)];
    inBetween = [WHO{1:length(results.timePoints),2}', fliplr(WHO{1:length(results.timePoints),6}')];
    fill(x2, inBetween*1000,[0.96 0.96 0.96],  'DisplayName', 'WHO 0-100% quartile');
    hold on
    inBetween2 = [WHO{1:length(results.timePoints),3}', fliplr(WHO{1:length(results.timePoints),5}')];
    fill(x2, inBetween2*1000,[0.91 0.91 0.92], 'DisplayName', 'WHO 25-75% quartile');
    hold on
    plot(results.timePoints,results.weight, 'LineWidth',4, 'DisplayName', 'infantWBM')
    hold on
    plot(results.timePoints,STIGMet{1:length(results.timePoints),1}, 'LineWidth',4, 'DisplayName', 'STIGMet')
    hold off

    set(gca,'box','off')
    xlabel('Age in days','FontSize',18)
    ylabel('Body weight (g)', 'FontSize',18,'FontName', 'Arial')
    title('Female infant growth prediction', 'FontSize',22,'FontName', 'Arial')
    legend('Location','northwest', 'FontSize', 13)
    
end
end