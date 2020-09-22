clear all; close all;
addpath(genpath('..\MATLAB Functions'))
%P0 Baseline
dir = './Data/P0_Snap25GC6s_MRS2500/*struct2.mat';
loadFileList(dir)
SGNstructs = loadCellStructs(dir)

%P0 base
dir = './Data/P0_Snap25GC6s_MRS2500/*struct2MRS2500.mat';
SGNstructsMRS2500 = loadCellStructs(dir);
%%
ROIsForGroup = 35; 
GroupBase_control = groupedActivity(SGNstructs,ROIsForGroup);
GroupBase_MRS = groupedActivity(SGNstructsMRS2500,ROIsForGroup);
%%

%%
for i=1:size(SGNstructs,2)
    if size([SGNstructs(i).events.isEmpty]) == size([SGNstructsMRS2500(i).events.isEmpty])
        masks = [SGNstructs(i).events.isEmpty] + [SGNstructsMRS2500(i).events.isEmpty];
        masks = masks >= 2;
        SGNstructs(i).meanFreq_ActiveD = mean([SGNstructs(i).events(~masks).frequency]);
        SGNstructsMRS2500(i).meanFreq_ActiveD = mean([SGNstructsMRS2500(i).events(~masks).frequency]);
    else
        i
    end
end

%% Testing for grouped activity

SGNstructsComb(1:2:size(SGNstructs,2)*2) = SGNstructs;
SGNstructsComb(2:2:size(SGNstructs,2)*2) = SGNstructsMRS2500;
for i = 1:50
     GroupBase = groupedActivity(SGNstructsComb,i);
     freqsBaseline(:,i) = [GroupBase(1:2:end).freq];
     freqsMRS(:,i) = [GroupBase(2:2:end).freq];
     roisBaseline(:,i) = [GroupBase(1:2:end).meanROIs];
     roisMRS(:,i) = [GroupBase(2:2:end).meanROIs];
 end

controlPrc = reshape([GroupBase_control.freqPrcActive],5,[])';
MRSPrc = reshape([GroupBase_MRS.freqPrcActive],5,[])';

figure; plot(mean(freqsBaseline,1)); hold on;
plot(mean(freqsMRS,1));

figure; plot(mean(roisBaseline,1)); hold on;
plot(mean(roisMRS,1));
%%
h = compare2P([SGNstructs.meanCorr_Active]', [SGNstructsMRS2500.meanCorr_Active]', {"Baseline","+MRS2500"}, 'Correlation coefficient', [3 2], [15 10]);
ylim([0 1]);
yticks([0:.25:1]);
xtickangle(45);
h1 = compare2P([SGNstructs.meanFreq_ActiveD]', [SGNstructsMRS2500.meanFreq_ActiveD]', {"Baseline","+MRS2500"}, 'Active ROI frequency', [3 2], [15 10]);
ylim([0 0.6]);
yticks(0:0.2:0.6);
xtickangle(45);
h2 = compare2P([SGNstructs.meanAmplitude_Active]'*100, [SGNstructsMRS2500.meanAmplitude_Active]'*100, {"Baseline","+MRS2500"}, 'Amplitude', [3 2], [15 10]);
ylim([0 300]);
%yticks(0:0.2:0.6);
xtickangle(45);
h3 = compare2P([SGNstructs.activeArea]', [SGNstructsMRS2500.activeArea]', {"Baseline","+MRS2500"}, 'Active area (%)', [3 2], [15 10]);
ylim([0 1]);
yticks([0:.25:1]);
xtickangle(45);
h4 = compare2P([SGNstructs.meanHW_Active]', [SGNstructsMRS2500.meanHW_Active]', {"Baseline","+MRS2500"}, 'Half-width (s)', [3 2], [15 10]);
ylim([0 12])
yticks([0:3:12]);
xtickangle(45);
h5 = compare2P([GroupBase_control.freq]', [GroupBase_MRS.freq]', {"Baseline","+MRS2500"}, 'Correlated events per min', [3 2], [15 10]);
ylim([0 3]);
yticks([0:1:3]);
xtickangle(45);
h6 = compare2P([GroupBase_control.meanROIs]', [GroupBase_MRS.meanROIs]', {"Baseline","+MRS2500"}, '# of ROIs per event', [3 2], [15 10]);
ylim([0 300]);
yticks([0:100:300]);
xtickangle(45);
handleTheSubplot({h5,h,h3,h1},[1 4])
figQuality(gcf,gca,[4.2 1.68]);
export_fig('.\EPS Panels\P0_MRS2500_quantification.eps');

%% statistics

%%

disp(['n = ' num2str(size(SGNstructs,2)/2)]);
[~,p, ci, stats] = ttest([GroupBase_control.freq]', [GroupBase_MRS.freq]');
disp(['Correlated events per minute p-value = ' num2str(p) '  t-stat: ' num2str(stats.tstat) '    tdf: ' num2str(stats.df)]);

[~,p,~,stats] = ttest([SGNstructs.meanCorr_Active]', [SGNstructsMRS2500.meanCorr_Active]');
disp(['Corr p-value = ' num2str(p) '  t-stat: ' num2str(stats.tstat) '    tdf: ' num2str(stats.df)]);

[~,p,~,stats] = ttest([SGNstructs.activeArea]', [SGNstructsMRS2500.activeArea]');
disp(['Active Roi area (%) p-value = ' num2str(p) '  t-stat: ' num2str(stats.tstat) '    tdf: ' num2str(stats.df)]);

[~,p,~,stats] = ttest([SGNstructs.meanFreq_ActiveD]', [SGNstructsMRS2500.meanFreq_ActiveD]');
disp(['Freq p-value = ' num2str(p) '  t-stat: ' num2str(stats.tstat) '    tdf: ' num2str(stats.df)]);

%not used
% [~,p] = ttest([conSGNstructs.meanAmplitude_Active]', [SGNstructsMRS2500.meanAmplitude_Active]');
% disp(['Amplitude p-value = ' num2str(p)]);
% 
% [~,p] = ttest([conSGNstructs.meanHW_Active]', [SGNstructsMRS2500.meanHW_Active]');
% disp(['HW p-value = ' num2str(p)]);


numTests = 6;
disp(['Bonferroni correction = 0.05 / ' num2str(numTests) ' = ' num2str(0.05/numTests)]);


%%
%%
% cellNum = 4
% rois = conSGNstructs(cellNum).rois;
% roisCNQX = cppSGNstructs(cellNum).rois;
% [t,n] = size(rois);
% 
% numToShow = 100;
% rng(1)
% temp = randperm(n)';
% temp = sort(temp(1:numToShow));
% isEmptyCon = [conSGNstructs(cellNum).events(temp).isEmpty];
% isEmptyCpp = [cppSGNstructs(cellNum).events(temp).isEmpty];
% offsets = -0.2*repmat(1:numToShow,t,1);
% 
% for i = 1:numToShow
%     rois(:,i) = smooth(rois(:,temp(i)),3);
%     roisCNQX(:,i) = smooth(roisCNQX(:,temp(i)),3);
% end
% 
% figure; plot(rois(1:275,isEmptyCon)+offsets(1:275,isEmptyCon),'Color',[0.7 0.7 0.7]); hold on;
% plot(rois(1:275,~isEmptyCon)+offsets(1:275,~isEmptyCon),'Color','k');
% ylim([-21 5]);
% xlim([0 275]);
% figQuality(gcf,gca,[2.4 1.6])
% export_fig('.\EPS Panels\baseline_CNQX_example.eps')
% % 
% figure; plot(roisCNQX(:,isEmptyCpp)+offsets(1:275,isEmptyCpp),'Color',[0.7 0.7 0.7]); hold on;
% plot(roisCNQX(:,~isEmptyCpp)+offsets(1:275,~isEmptyCpp),'Color','k');
% ylim([-21 5]);
% xlim([0 275]);
% figQuality(gcf,gca,[2.4 1.6])
% export_fig('.\EPS Panels\CNQX_example.eps')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%example images
%cochlea_7_base_1 and base_1_MRS2500 
%baseline 601-875, MRS 301-575
%median filt 3d (1x1x1) then normalize contrast (18 to 3400), in
%data/Time_color_Processed
%%Example traces
cellNum = 7
rois = SGNstructs(cellNum).rois;
roisMRS = SGNstructsMRS2500(cellNum).rois;
[t,n] = size(rois)

numToShow = 100;
rng(3)
temp = randperm(n)'
temp = sort(temp(1:numToShow));
isEmptyCon = [SGNstructs(cellNum).events(temp).isEmpty];
isEmptyMRS = [SGNstructsMRS2500(cellNum).events(temp).isEmpty];
offsets = -0.2*repmat(1:numToShow,t,1);

for i = 1:numToShow
    rois(:,i) = smooth(rois(:,temp(i)),3);
    roisMRS(:,i) = smooth(roisMRS(:,temp(i)),3);
end

figure; plot(rois(601:875,isEmptyCon)+offsets(601:875,isEmptyCon),'Color',[0.7 0.7 0.7]); hold on;
plot(rois(601:875,~isEmptyCon)+offsets(601:875,~isEmptyCon),'Color','k');
ylim([-21 7]);
xlim([0 275]);
figQuality(gcf,gca,[2.4 1.6]);
export_fig('.\EPS Panels\baseline_example.eps')

figure; plot(roisMRS(301:575,isEmptyMRS)+offsets(601:875,isEmptyMRS),'Color',[0.7 0.7 0.7]); hold on;
plot(roisMRS(301:575,~isEmptyMRS)+offsets(601:875,~isEmptyMRS),'Color','k');
ylim([-21 7]);
xlim([0 275]);
figQuality(gcf,gca,[2.4 1.6]);
export_fig('.\EPS Panels\MRS_example.eps')
%% E16.5

clear all; close all;
%E16.5 Baseline
dir = './Data/E16.5_MRS2500_1/*struct2.mat';
SGNstructs = loadCellStructs(dir)
dir = './Data/E16.5_MRS2500_2/*struct2.mat';
SGNstructs = [SGNstructs loadCellStructs(dir)];

%E16.5 base
dir = './Data/E16.5_MRS2500_1/*struct2MRS2500.mat';
SGNstructsMRS2500 = loadCellStructs(dir);
dir = './Data/E16.5_MRS2500_2/*struct2MRS2500.mat';
SGNstructsMRS2500 = [SGNstructsMRS2500 loadCellStructs(dir)];
%%
ROIsForGroup = 35; 
GroupBase_control = groupedActivity(SGNstructs,ROIsForGroup);
GroupBase_MRS = groupedActivity(SGNstructsMRS2500,ROIsForGroup);

%% Testing for grouped activity

SGNstructsComb(1:2:size(SGNstructs,2)*2) = SGNstructs;
SGNstructsComb(2:2:size(SGNstructs,2)*2) = SGNstructsMRS2500;
for i = 1:50
     GroupBase = groupedActivity(SGNstructsComb,i);
     freqsBaseline(:,i) = [GroupBase(1:2:end).freq];
     freqsMRS(:,i) = [GroupBase(2:2:end).freq];
     roisBaseline(:,i) = [GroupBase(1:2:end).meanROIs];
     roisMRS(:,i) = [GroupBase(2:2:end).meanROIs];
 end

controlPrc = reshape([GroupBase_control.freqPrcActive],5,[])';
MRSPrc = reshape([GroupBase_MRS.freqPrcActive],5,[])';

figure; plot(mean(freqsBaseline,1)); hold on;
plot(mean(freqsMRS,1));

figure; plot(mean(roisBaseline,1)); hold on;
plot(mean(roisMRS,1));
%%
h = compare2P([SGNstructs.meanCorr_Active]', [SGNstructsMRS2500.meanCorr_Active]', {"Baseline","+MRS2500"}, 'Corr', [3 2], [10 15]);
ylim([0 1]);
yticks([0:.25:1]);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h1 = compare2P([SGNstructs.meanFreq_Active]', [SGNstructsMRS2500.meanFreq_Active]', {"Baseline","+MRS2500"}, 'Freq', [3 2], [10 15]);
ylim([0 0.6]);
yticks(0:0.2:0.6);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h2 = compare2P([SGNstructs.meanAmplitude_Active]'*100, [SGNstructsMRS2500.meanAmplitude_Active]'*100, {"Baseline","+MRS2500"}, 'Amplitude', [3 2], [10 15]);
ylim([0 300]);
%yticks(0:0.2:0.6);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h3 = compare2P([SGNstructs.activeArea]', [SGNstructsMRS2500.activeArea]', {"Baseline","+MRS2500"}, 'Active area (%)', [3 2], [10 15]);
ylim([0 1]);
yticks([0:.25:1]);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h4 = compare2P([SGNstructs.meanHW_Active]', [SGNstructsMRS2500.meanHW_Active]', {"Baseline","+MRS2500"}, 'Half-width (s)', [3 2], [10 15]);
ylim([0 12])
yticks([0:3:12]);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h5 = compare2P([GroupBase_control.freq]', [GroupBase_MRS.freq]', {"Baseline","+MRS2500"}, 'Correlated events per minute', [3 2], [10 15]);
ylim([0 2]);
yticks([0:1:2]);
xtickangle(45);
h6 = compare2C([GroupBase_control.meanROIs]', [GroupBase_MRS.meanROIs]', {"Baseline","+MRS2500"}, '# of ROIs per event', [3 2], [10 15]);
ylim([0 300]);
yticks([0:100:300]);
xtickangle(45);
handleTheSubplot({h3,h1,h2,h4,h,h5,h6},[1 7]);
figQuality(gcf,gca,[9.25 1.5]);
export_fig('.\EPS Panels\P0_MRS2500_quantification.eps');

%% statistics

[~,p] = ttest([SGNstructs.meanCorr_Active]', [SGNstructsMRS2500.meanCorr_Active]');

disp(['n = ' num2str(size(SGNstructs,2))]);
disp(['Corr p-value = ' num2str(p)]);
[~,p] = ttest([SGNstructs.meanFreq_Active]', [SGNstructsMRS2500.meanFreq_Active]');
disp(['Freq p-value = ' num2str(p)]);
[~,p] = ttest([SGNstructs.meanAmplitude_Active]', [SGNstructsMRS2500.meanAmplitude_Active]');
disp(['Amplitude p-value = ' num2str(p)]);
[~,p] = ttest([SGNstructs.activeArea]', [SGNstructsMRS2500.activeArea]');
disp(['Active Roi area (%) p-value = ' num2str(p)]);
[~,p] = ttest([SGNstructs.meanHW_Active]', [SGNstructsMRS2500.meanHW_Active]');
disp(['HW p-value = ' num2str(p)]);
[~,p] = ttest([GroupBase_control.freq]', [GroupBase_MRS.freq]');
disp(['Correlated event p-value = ' num2str(p)]);
[~,p] = ttest([GroupBase_control.meanROIs]', [GroupBase_MRS.meanROIs]');
disp(['Correlated ROIs p-value = ' num2str(p)]);

numTests = 6;
disp(['Bonferroni correction = 0.05 / ' num2str(numTests) ' = ' num2str(0.05/numTests)]);