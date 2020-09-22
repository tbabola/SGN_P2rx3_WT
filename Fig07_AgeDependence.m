%%
clear all; close all;
addpath(genpath('..\MATLAB Functions'))

%%
img = loadTif('.\Data\P0 Base Example\STD_cochlea_1_base_2.tif',32);
img = img(:,:,1);
imshow(img);
greencmap = zeros(size(gfb));
greencmap(:,2) = [0:1:255]/255;
colormap(greencmap)
caxis([0 1400])

%%

SGNStruct = load('.\Data\P0 Base Example\cochlea_1_base_2_SGNstruct2.mat');
SGNStruct = SGNStruct.SGNstruct;

widthImg = size(img,2);
heightImg = size(img,1);
sizeSq = 10;

indices = getGrid(widthImg,heightImg,sizeSq);
positiveIndices = getPositiveGrid(indices,SGNStruct.SGNmask.Position);
positiveIndices(:,9:10) = positiveIndices(:,1:2);
for i=1:size(positiveIndices,1)
    hold on;
    plot(positiveIndices(i,1:2:end),positiveIndices(i,2:2:end),'Color','g');
end

export_fig('.\EPS Panels\example_SGN_grid.eps');

%%
rois = SGNStruct.rois;
n = size(rois,2);
t = 300;
numToShow = 100;
rng(16);
temp = randperm(n)'
temp = sort(temp(1:numToShow));
isEmpty = [SGNStruct.events(temp).isEmpty];
offsets = -0.2*repmat(1:numToShow,t,1);
plotrois = rois(301:600,temp);
for i = 1:numToShow
    plotrois(:,i) = smooth(plotrois(:,i),3);
end
figure; plot(plotrois(:,isEmpty)+offsets(:,isEmpty),'Color',[0.7 0.7 0.7]); hold on;
plot(plotrois(:,~isEmpty)+offsets(:,~isEmpty),'Color','k'); 
figQuality(gcf,gca,[3.1 1.9]);
ylim([-20 10])
export_fig('.\EPS Panels\example_ROIs.eps');

%%
%Note: time-color figures were generated in FIJI. Subtract off mean of
%movie, median filter 3d (2 2 2), set limits to 0 to 1400, run time-color histogram

%%
%E16.5 Apex
SGNStruct = load('C:\Users\Bergles Lab\Desktop\SGN spontaneous activity for Travis\SGN spontaneous activity for Travis\E16.5 snap25-gcamp6s_CNQX\apex\cochlea_4_apex_1_SGNstruct2.mat');
SGNStruct = SGNStruct.SGNstruct;
rois = SGNStruct.rois;
n = size(rois,2);
t = 300;
numToShow = 100;
rng(2);
temp = randperm(n)'
temp = sort(temp(1:numToShow));
isEmpty = [SGNStruct.events(temp).isEmpty];
offsets = -0.2*repmat(1:numToShow,t,1);
plotrois = rois(1:300,temp);
for i = 1:numToShow
    plotrois(:,i) = smooth(plotrois(:,i),3);
end

figure; plot(plotrois(:,isEmpty)+offsets(:,isEmpty),'Color',[0.7 0.7 0.7]); hold on;
plot(plotrois(:,~isEmpty)+offsets(:,~isEmpty),'Color','k'); 
figQuality(gcf,gca,[2.6 1.9]);
ylim([-20 5])
export_fig('.\EPS Panels\E165_ROIs.eps');

%%
SGNStruct = load('C:\Users\Bergles Lab\Desktop\SGN spontaneous activity for Travis\SGN spontaneous activity for Travis\E16.5 snap25-gcamp6s_CNQX\base\cochlea_3_base_1_SGNstruct2.mat');
SGNStruct = SGNStruct.SGNstruct;
rois = SGNStruct.rois;
n = size(rois,2);
t = 300;
numToShow = 100;
rng(11);
temp = randperm(n)'
temp = sort(temp(1:numToShow));
isEmpty = [SGNStruct.events(temp).isEmpty];
offsets = -0.2*repmat(1:numToShow,t,1);
plotrois = rois(1:300,temp);
for i = 1:numToShow
    plotrois(:,i) = smooth(plotrois(:,i),3);
end

figure; plot(plotrois(:,isEmpty)+offsets(:,isEmpty),'Color',[0.7 0.7 0.7]); hold on;
plot(plotrois(:,~isEmpty)+offsets(:,~isEmpty),'Color','k'); 
figQuality(gcf,gca,[2.6 1.9]);
ylim([-20 5])
export_fig('.\EPS Panels\E165_base_ROIs.eps');

%%
%%
SGNStruct = load('C:\Users\Bergles Lab\Desktop\SGN spontaneous activity for Travis\SGN spontaneous activity for Travis\P0 WT\apex\cochlea_4_apex_1_SGNstruct2.mat');
SGNStruct = SGNStruct.SGNstruct;
rois = SGNStruct.rois;
n = size(rois,2);
t = 300;
numToShow = 100;
rng(11);
temp = randperm(n)'
temp = sort(temp(1:numToShow));
isEmpty = [SGNStruct.events(temp).isEmpty];
offsets = -0.2*repmat(1:numToShow,t,1);
plotrois = rois(1:300,temp);
for i = 1:numToShow
    plotrois(:,i) = smooth(plotrois(:,i),3);
end

figure; plot(plotrois(:,isEmpty)+offsets(:,isEmpty),'Color',[0.7 0.7 0.7]); hold on;
plot(plotrois(:,~isEmpty)+offsets(:,~isEmpty),'Color','k'); 
figQuality(gcf,gca,[2.6 1.9]);
ylim([-20 5])
export_fig('.\EPS Panels\P0_apex_ROIs.eps');

%%
%%
SGNStruct = load('C:\Users\Bergles Lab\Desktop\SGN spontaneous activity for Travis\SGN spontaneous activity for Travis\P0 WT\base\cochlea_4_base_2_SGNstruct2.mat');
SGNStruct = SGNStruct.SGNstruct;
rois = SGNStruct.rois;
n = size(rois,2);
t = 300;
numToShow = 100;
rng(1);
temp = randperm(n)'
temp = sort(temp(1:numToShow));
isEmpty = [SGNStruct.events(temp).isEmpty];
offsets = -0.2*repmat(1:numToShow,t,1);
plotrois = rois(1:300,temp);
for i = 1:numToShow
    plotrois(:,i) = smooth(plotrois(:,i),3);
end

figure; plot(plotrois(:,isEmpty)+offsets(:,isEmpty),'Color',[0.7 0.7 0.7]); hold on;
plot(plotrois(:,~isEmpty)+offsets(:,~isEmpty),'Color','k'); 
figQuality(gcf,gca,[2.6 1.9]);
ylim([-20 5])
export_fig('.\EPS Panels\P0_base_ROIs.eps');

%%
%E16.5 to P0 transition
basedir = 'C:\Users\Bergles Lab\Desktop\SGN spontaneous activity for Travis\SGN spontaneous activity for Travis\';

%dir = './SGN spontaneous activity for Travis/E16.5 SGN spontaneous activity_snap25-GCaMP6s/apex/*struct2.mat';
EapexSGNstructs = loadCellStructs([basedir 'E16.5_WT*/apex/*SGNstruct2.mat']);
EbaseSGNstructs = loadCellStructs([basedir 'E16.5_WT_s*/base/*SGNstruct2.mat']);
EbaseSGNstructs = [EbaseSGNstructs loadCellStructs([basedir 'E16.5_WT_M*/Old Analysis Backup/*SGNstruct2.mat'])];
PbaseSGNstructs = loadCellStructs([basedir 'P0 WT*/base/*SGNstruct2.mat']);
PapexSGNstructs = loadCellStructs([basedir 'P0 WT*/apex/*SGNstruct2.mat']);

%%
h = compare4_jitter([EapexSGNstructs.meanCorr_Active]', [EbaseSGNstructs.meanCorr_Active]', [PapexSGNstructs.meanCorr_Active]',[PbaseSGNstructs.meanCorr_Active]', {"Apex","Base","Apex","Base"}, 'Corr', [3 2], [5 10]);
ylim([0 1]);
yticks([0:.25:1]);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h1 = compare4_jitter([EapexSGNstructs.meanFreq_Active]', [EbaseSGNstructs.meanFreq_Active]', [PapexSGNstructs.meanFreq_Active]',[PbaseSGNstructs.meanFreq_Active]', {"Apex","Base","Apex","Base"}, 'Freq', [3 2], [5 10]);
ylim([0 0.6]);
yticks(0:0.2:0.6);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h2 = compare4_jitter([EapexSGNstructs.meanAmplitude_Active]', [EbaseSGNstructs.meanAmplitude_Active]', [PapexSGNstructs.meanAmplitude_Active]',[PbaseSGNstructs.meanAmplitude_Active]', {"Apex","Base","Apex","Base"}, 'Amplitude', [3 2], [5 10]);
ylim([0 3]);
%yticks(0:0.2:0.6);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h3 = compare4_jitter([EapexSGNstructs.activeArea]', [EbaseSGNstructs.activeArea]', [PapexSGNstructs.activeArea]',[PbaseSGNstructs.activeArea]', {"Apex","Base","Apex","Base"}, 'Active area (%)', [3 2], [5 10]);
ylim([0 1]);
yticks([0:.25:1]);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);
h4 = compare4_jitter([EapexSGNstructs.meanHW_Active]', [EbaseSGNstructs.meanHW_Active]', [PapexSGNstructs.meanHW_Active]',[PbaseSGNstructs.meanHW_Active]', {"Apex","Base","Apex","Base"}, 'Half-width (s)', [3 2], [5 10]);
ylim([0 12])
yticks([0:3:12]);
xtickangle(45);
figQuality(gcf,gca,[2 1.5]);

%%
ROIsForGroup = 35; 
EGroupBase = groupedActivity(EbaseSGNstructs,ROIsForGroup);
EGroupApex = groupedActivity(EapexSGNstructs,ROIsForGroup);
PGroupBase = groupedActivity(PbaseSGNstructs,ROIsForGroup);
PGroupApex = groupedActivity(PapexSGNstructs,ROIsForGroup);

h5 = compare4_jitter([EGroupApex.freq]', [EGroupBase.freq]', [PGroupApex.freq]',[PGroupBase.freq]', {"Apex","Base","Apex","Base"}, 'Correlated events per minute', [3 2], [5 10]);
ylim([0 4]);
yticks([0:1:4]);
xtickangle(45);
h6 = compare4_jitter([EGroupApex.meanROIs]', [EGroupBase.meanROIs]', [PGroupApex.meanROIs]',[PGroupBase.meanROIs]', {"Apex","Base","Apex","Base"}, 'Active ROIs per event', [3 2], [5 10]);
ylim([0 250]);
yticks([0:50:250]);
xtickangle(45);
%handleTheSubplot({h,h1,h2,h3,h4,h5,h6},[1 7]);
%handleTheSubplot({h3,h,h5,h6,h1,h2,h4},[1 7]);
handleTheSubplot({h3,h,h5,h1,h2,h4},[1 6]);
figQuality(gcf,gca,[9 1.5]);
export_fig('.\EPS Panels\transient_quantification.eps')


%%
%STATISTICS
group = [ones(size(EGroupApex,2),1); 2*ones(size(EGroupBase,2),1); 3*ones(size(PGroupApex,2),1); 4*ones(size(PGroupBase,2),1)];

activeArea = [[EapexSGNstructs.activeArea]'; [EbaseSGNstructs.activeArea]'; [PapexSGNstructs.activeArea]';[PbaseSGNstructs.activeArea]'];
[p,anovatab1,stats]=anova1(activeArea,group);
c = multcompare(stats);

corrSGN = [[EapexSGNstructs.meanCorr_Active]'; [EbaseSGNstructs.meanCorr_Active]'; [PapexSGNstructs.meanCorr_Active]';[PbaseSGNstructs.meanCorr_Active]'];
[p1,anovatab2,stats]=anova1(corrSGN,group);
c1 = multcompare(stats);

coordFreq = [[EGroupApex.freq]'; [EGroupBase.freq]'; [PGroupApex.freq]'; [PGroupBase.freq]'];
[p2,anovatab3,stats]=anova1(coordFreq,group);
c2 = multcompare(stats);

freqSGN = [[EapexSGNstructs.meanFreq_Active]'; [EbaseSGNstructs.meanFreq_Active]'; [PapexSGNstructs.meanFreq_Active]';[PbaseSGNstructs.meanFreq_Active]'];
[p3,anovatab4,stats]=anova1(freqSGN,group);
c3 = multcompare(stats);

ampSGN = [[EapexSGNstructs.meanAmplitude_Active]'; [EbaseSGNstructs.meanAmplitude_Active]'; [PapexSGNstructs.meanAmplitude_Active]';[PbaseSGNstructs.meanAmplitude_Active]'];
[p4,anovatab5,stats]=anova1(ampSGN,group);
c4 = multcompare(stats);

hwSGN = [[EapexSGNstructs.meanHW_Active]'; [EbaseSGNstructs.meanHW_Active]'; [PapexSGNstructs.meanHW_Active]';[PbaseSGNstructs.meanHW_Active]'];
[p5,anovatab6,stats]=anova1(hwSGN,group);
c5 = multcompare(stats);

coordROIs = [[EGroupApex.meanROIs]'; [EGroupBase.meanROIs]'; [PGroupApex.meanROIs]'; [PGroupBase.meanROIs]'];
[p6,anovatab7,stats]=anova1(coordROIs,group);
c6 = multcompare(stats);