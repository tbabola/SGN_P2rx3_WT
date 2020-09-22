%%
% load baseline files and refind peaks for 10 minute period (300 to 900s)
files = loadFileList('.\Data\E16.5_MRS2500_2\*SGNstruct2.mat');

for i = 1:size(files,1)
    h = load(files{i});
    SGNstructO = h.SGNstruct;
    SGNstruct = analyzeROIs(SGNstructO.rois(300:end,:));
    SGNstruct.rois = SGNstructO.rois;
    save(files{i},'SGNstruct');
end

%% load MRS files and refind peaks for 10 minute period (300 to 900s)
filesMRS = loadFileList('.\Data\E16.5_MRS2500_2\*SGNstruct2MRS2500.mat');
for i = 1:size(filesMRS,1)
    h = load(filesMRS{i});
    p = load(files{i});
    filesMRS{i}
    files{i}
    SGNstructO = h.SGNstruct;
    SGNstructbaseline = p.SGNstruct;
    SGNstruct = analyzeROIs(SGNstructO.rois(300:end,:), SGNstructbaseline.stds);
    SGNstruct.rois = SGNstructO.rois;
    save(filesMRS{i},'SGNstruct');
end
%%

%%


function SGNstruct = analyzeROIs(rois,stds)
    if nargin == 1
        rois_std = std(rois,[],1);
    else
        rois_std = stds;
    end
    rois_std(rois_std >  0.1) =  0.1;
    rois_std(rois_std <  prctile(rois_std,5)) =  prctile(rois_std,5);
    [t,n] = size(rois);
    rois_mean = mean(rois,1);
    events = struct();
    timeMinutes = size(rois,1)/60;
    SGNstruct.corrROIs = corr(rois);
    SGNstruct.corrshuffROIs = corr(shuffleSignals(rois));
    warning('off','signal:findpeaks:largeMinPeakHeight'); %suppress warning that states there are no peaks detected
    binarywhole = [];
    parfor i=1:n
        %[pks,locs,w] = findpeaks(rois_adj(:,i),1:t);
        [pks,locs,w] = findpeaks(rois(:,i),1:t,'MinPeakHeight',rois_mean(i)+5*rois_std(i));
        events(i).pks = pks;
        events(i).locs = locs;
        temp = zeros(1,t);
        temp(locs) = 1;
        temp = conv(temp,[1 1 1],'same');
        events(i).binary = temp;
        events(i).w = w;
        events(i).isEmpty = isempty(pks);
        events(i).meanAmp = mean(pks);
        events(i).frequency = size(locs,2)/timeMinutes;
        events(i).corr80 = prctile(SGNstruct.corrROIs(:,i),80,1);
        events(i).corrshuf80 = prctile(SGNstruct.corrshuffROIs(:,i),80,1);
        binarywhole(i,:) = temp;
    end
    
    SGNstruct.means = rois_mean;
    SGNstruct.stds = rois_std;
    SGNstruct.activeArea = 1 - sum([events.isEmpty])/size(events,2);
    SGNstruct.totalROIs = size(events,2);
    SGNstruct.events = events;
    SGNstruct.meanAmplitude_Active = mean([events(~[events.isEmpty]).meanAmp]);
    SGNstruct.meanFreq_Active = mean([events(~[events.isEmpty]).frequency]);
    SGNstruct.meanHW_Active = mean([events(~[events.isEmpty]).w]);
    SGNstruct.meanFreq_AllROIs = mean([events.frequency]);
    SGNstruct.meanCorr_Active = mean([events(~[events.isEmpty]).corr80]);
    SGNstruct.meanshufCorr_Active = mean([events(~[events.isEmpty]).corrshuf80]);
end