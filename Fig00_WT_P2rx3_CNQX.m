addpath(genpath("MATLAB Functions"))
%% P0 base WT CNQX/CPP Experiments
clear all; close all;

dir = './Data/P0 WT*/base/*struct2*.mat';
files = findPairs(dir,'cpp');
SGNstructs = loadCellStructs(files);
h = analyzeAndPlot(SGNstructs);
export_fig('.\EPS Panels\P0_WT_CNQX_CPP.eps');

%% P0 base P2rx3 CNQX/CPP Experiments
clear all; close all;
dir = './Data/P0 P2rx3 KO*/base/*struct2*.mat';
files = findPairs(dir,'cpp');
SGNstructs = loadCellStructs(files);
h = analyzeAndPlot(SGNstructs);
export_fig('.\EPS Panels\P0_P2rx3KO_CNQX_CPP.eps');

%% P0 WT versus P2rx3 KO
clear all; close all;
dir = './Data/P0 WT*/base/*struct2.mat';
SGNstructs = loadCellStructs(dir);

dir = './Data/P0 P2rx3 KO*/base/*struct2.mat';
SGNstructsKO = loadCellStructs(dir);

h = analyzeAndPlotWTvKO(SGNstructs, SGNstructsKO);
export_fig('.\EPS Panels\P0_WTvKO.eps');

%% E16.5 base WT CNQX/CPP Experiments
clear all; close all;

dir = './Data/E16.5_WT*_CNQX_CPP*/base/*struct2*.mat';
files = findPairs(dir,'cpp');
SGNstructs = loadCellStructs(files);
h = analyzeAndPlot(SGNstructs);
export_fig('.\EPS Panels\E165_WT_CNQX_CPP.eps');

%% E16.5 base P2rx3 CNQX/CPP Experiments
clear all; close all;

dir = './Data/E16.5 P2rx3 KO*/base/*struct2*.mat';
files = findPairs(dir,'cpp');
SGNstructs = loadCellStructs(files);
h = analyzeAndPlot(SGNstructs);
export_fig('.\EPS Panels\E165_P2rx3KO_CNQX_CPP.eps');

%% E16.5 WT versus P2rx3 KO
clear all; close all;
dir = './Data/E16.5_WT*/base/*struct2.mat';
fileList = loadFileList(dir);
for i=1:size(fileList,1)
    s = load(fileList{i});
    fns = fieldnames(s);
    temp = s.(fns{1});
    temp.rotateDegrees = 0;
    temp.SGNmask = [];
    tempCells(i)= temp;
end
SGNstructs = tempCells;

dir = './Data/E16.5 P2rx3 KO*/base/*struct2.mat';
SGNstructsKO = loadCellStructs(dir);

h = analyzeAndPlotWTvKO(SGNstructs, SGNstructsKO);
export_fig('.\EPS Panels\E165_WTvKO.eps');
%%
cellNum = 2;
rois = conSGNstructs(cellNum).rois;
roisCNQX = cppSGNstructs(cellNum).rois;
[t,n] = size(rois);

numToShow = 100;
rng(1)
temp = randperm(n)'
temp = sort(temp(1:numToShow));
for i = 1:numToShow
    rois(:,i) = smooth(rois(:,i),3);
    roisCNQX(:,i) = smooth(roisCNQX(:,i),3);
end

figure; plot(rois(1:300,temp) - 0.2*repmat(1:numToShow,300,1),'Color','k');
ylim([-20 5]);
figQuality(gcf,gca,[2.6 1.5])
export_fig('.\EPS Panels\E165_baseline_CNQX_example.eps')

figure; plot(roisCNQX(:,temp) - 0.2*repmat(1:numToShow,size(roisCNQX,1),1),'Color','k');
ylim([-20 5])
figQuality(gcf,gca,[2.6 1.5])
export_fig('.\EPS Panels\E165_CNQX_example.eps')


%%
function h = analyzeAndPlot(SGNstructs)
    %go through the structs and find only cells that are active in either
    %the baseline or drug condition for meanFreq calculation
    for i=1:2:size(SGNstructs,2)
        if size([SGNstructs(i).events.isEmpty]) == size([SGNstructs(i+1).events.isEmpty])
            masks = [SGNstructs(i).events.isEmpty] + [SGNstructs(i+1).events.isEmpty];
            masks = masks >= 2;
            SGNstructs(i).meanFreq_ActiveD = mean([SGNstructs(i).events(~masks).frequency]);
            SGNstructs(i+1).meanFreq_ActiveD = mean([SGNstructs(i+1).events(~masks).frequency]);
        else
            i
        end
    end
    
    %set up and analyze group activity
    ROIsForGroup = 35;
    GroupBase = groupedActivity(SGNstructs,ROIsForGroup);
    
    conSGNstructs = SGNstructs(1:2:end);
    cppSGNstructs = SGNstructs(2:2:end);
    CNQXcolor = [102 153 255]/255;
    h = compare2P([conSGNstructs.meanCorr_Active]', [cppSGNstructs.meanCorr_Active]', {"Baseline","+CNQX/CPP"}, 'Correlation coefficient', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 1]);
    yticks(0:0.25:1);
    h1 =compare2P([conSGNstructs.meanFreq_ActiveD]', [cppSGNstructs.meanFreq_ActiveD]', {"Baseline","+CNQX/CPP"}, 'Active ROI frequency', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 0.5]);
    yticks(0:0.1:.5);
    h2 = compare2P([conSGNstructs.meanAmplitude_Active]', [cppSGNstructs.meanAmplitude_Active]', {"Baseline","+CNQX/CPP"}, 'Amp', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 3]);
    yticks(0:1:3);
    h3 = compare2P([conSGNstructs.activeArea]', [cppSGNstructs.activeArea]', {"Baseline","+CNQX/CPP"}, 'Active area (%)', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 1]);
    yticks(0:.25:1);
    h4 = compare2P([conSGNstructs.meanHW_Active]', [cppSGNstructs.meanHW_Active]', {"Baseline","+CNQX/CPP"}, 'half-width (s)', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 10]);
    yticks(0:2.5:10);
    h5 = compare2P([GroupBase(1:2:end).freq]', [GroupBase(2:2:end).freq]', {"Baseline","+CNQX/CPP"}, 'Correlated events per min', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 3]);
    yticks(0:1:3);
    h = handleTheSubplot({h5,h1,h2,h3,h4,h},[1 6]);
    %h = handleTheSubplot({h5,h,h3,h1},[1 4])
    figQuality(gcf,gca,[6 1.75]);
    
    [~,p] = ttest([conSGNstructs.meanCorr_Active]', [cppSGNstructs.meanCorr_Active]');
    disp(['n = ' num2str(size(SGNstructs,2)/2)]);
    disp(['Corr p-value = ' num2str(p)]);
    [~,p] = ttest([conSGNstructs.meanFreq_ActiveD]', [cppSGNstructs.meanFreq_ActiveD]');
    disp(['Freq p-value = ' num2str(p)]);
    [~,p] = ttest([conSGNstructs.meanAmplitude_Active]', [cppSGNstructs.meanAmplitude_Active]');
    disp(['Amplitude p-value = ' num2str(p)]);
    [~,p] = ttest([conSGNstructs.activeArea]', [cppSGNstructs.activeArea]');
    disp(['Active Roi area (%) p-value = ' num2str(p)]);
    [~,p] = ttest([conSGNstructs.meanHW_Active]', [cppSGNstructs.meanHW_Active]');
    disp(['HW p-value = ' num2str(p)]);
    [~,p] = ttest([GroupBase(1:2:end).freq]', [GroupBase(2:2:end).freq]');
    disp(['Correlated event p-value = ' num2str(p)]);

    numTests = 6;
    disp(['Bonferroni correction = 0.05 / ' num2str(numTests) ' = ' num2str(0.05/numTests)]);
end

function h = analyzeAndPlotWTvKO(SGNstructs, SGNstructsKO)

    %set up and analyze group activity
    ROIsForGroup = 35;
    GroupBase = groupedActivity(SGNstructs,ROIsForGroup);
    GroupBaseKO = groupedActivity(SGNstructsKO,ROIsForGroup);
    
    CNQXcolor = [102 153 255]/255;
    h = compare2([SGNstructs.meanCorr_Active]', [SGNstructsKO.meanCorr_Active]', {"WT","P2rx3 KO"}, 'Correlation coefficient', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 1]);
    yticks(0:0.25:1);
    h1 =compare2([SGNstructs.meanFreq_Active]', [SGNstructsKO.meanFreq_Active]', {"WT","P2rx3 KO"}, 'Active ROI frequency', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    %ylim([0 0.75]);
    %yticks(0:0.25:.75);
    h2 = compare2([SGNstructs.meanAmplitude_Active]', [SGNstructsKO.meanAmplitude_Active]', {"WT","P2rx3 KO"}, 'Amp', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 3]);
    yticks(0:1:3);
    h3 = compare2([SGNstructs.activeArea]', [SGNstructsKO.activeArea]', {"WT","P2rx3 KO"}, 'Active area (%)', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 1]);
    yticks(0:.25:1);
    h4 = compare2([SGNstructs.meanHW_Active]', [SGNstructsKO.meanHW_Active]', {"WT","P2rx3 KO"}, 'half-width (s)', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 10]);
    yticks(0:2.5:10);
    h5 = compare2([GroupBase(1:2:end).freq]', [GroupBase(2:2:end).freq]', {"WT","P2rx3 KO"}, 'Correlated events per min', [3 2], [15 10],'k',CNQXcolor);
    xtickangle(45);
    ylim([0 3]);
    yticks(0:1:3);
    h = handleTheSubplot({h5,h1,h2,h3,h4,h},[1 6]);
    %h = handleTheSubplot({h5,h,h3,h1},[1 4])
    figQuality(gcf,gca,[6 1.75]);
    
    [~,p] = ttest2([SGNstructs.meanCorr_Active]', [SGNstructsKO.meanCorr_Active]');
    disp(['n = ' num2str(size(SGNstructs,2)) ' WT and n = ' num2str(size(SGNstructsKO,2)) ' P2rx3 KO.']);
    disp(['Corr p-value = ' num2str(p)]);
    [~,p] = ttest2([SGNstructs.meanFreq_Active]', [SGNstructsKO.meanFreq_Active]');
    disp(['Freq p-value = ' num2str(p)]);
    [~,p] = ttest2([SGNstructs.meanAmplitude_Active]', [SGNstructsKO.meanAmplitude_Active]');
    disp(['Amplitude p-value = ' num2str(p)]);
    [~,p] = ttest2([SGNstructs.activeArea]', [SGNstructsKO.activeArea]');
    disp(['Active Roi area (%) p-value = ' num2str(p)]);
    [~,p] = ttest2([SGNstructs.meanHW_Active]', [SGNstructsKO.meanHW_Active]');
    disp(['HW p-value = ' num2str(p)]);
    [~,p] = ttest2([GroupBase(1:2:end).freq]', [GroupBase(2:2:end).freq]');
    disp(['Correlated event p-value = ' num2str(p)]);

    numTests = 6;
    disp(['Bonferroni correction = 0.05 / ' num2str(numTests) ' = ' num2str(0.05/numTests)]);
end