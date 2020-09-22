%%
%load image
clear all; close all;
[fn, dname] = uigetfile();
[dname fn]
baselineexists = 0;
drug = 'cpp'
if contains(fn,drug)
    baselineexists = 1;
    [tempfn, tempdname] = uigetfile(dname);
    blStruct = loadCellStructs([tempdname tempfn]);
end

SGNstruct = struct();

imgdata = bfopen([dname fn])
imgdata = imgdata{1,1};
img = []; index = 1;
for i = 1:2:size(imgdata,1)
    img(:,:,index) = imgdata{i,1};
    index = index + 1;
end
[fp,name,~] = fileparts([dname fn]);
%%
if ~baselineexists
    figure; imshow(mean(img,3)/mean(max(mean(img,3))));
    x = input('Rotate? (degrees)')
    SGNstruct.rotateDegrees = x;
    img = imrotate(img,x);
else
    img = imrotate(img,blStruct.rotateDegrees);
    SGNstruct.rotateDegrees = blStruct.rotateDegrees;
end
imshow(mean(img,3)/mean(max(mean(img,3)))); truesize;
%%
widthImg = size(img,2);
heightImg = size(img,1);
sizeSq = 10;

indices = getGrid(widthImg,heightImg,sizeSq);

if ~baselineexists
    roi = drawpolygon;
    SGNstruct.SGNmask = roi; %%save outline to struct
    roiIndices = roi.Position;
else
    roi = blStruct.SGNmask;
    SGNstruct.SGNmask = roi;
    roiIndices = roi.Position;
end
positiveIndices = getPositiveGrid(indices,roiIndices);

%for i=1:size(indices,1)
  %  hold on;
  %  plot(indices(i,1:2:end),indices(i,2:2:end),'Color','k');
%end

for i=1:size(positiveIndices,1)
    hold on;
    plot(positiveIndices(i,1:2:end),positiveIndices(i,2:2:end),'Color','g');
end
if ~baselineexists
    saveas(gcf,[fp '\' name '_gridImage.bmp'])
end
%%
%bleachCorrect
img = bleachCorrect(img,1);
rois = normalizeGridImg(img,10,positiveIndices);


[t,n] = size(rois);
SGNstruct.rois = rois;

figure; imagesc(corr(rois))
SGNstruct.corrROIs = corr(rois);
SGNstruct.corrshuffROIs = corr(shuffleSignals(rois));

numToShow = 200;
temp = randperm(n)'
temp = sort(temp(1:numToShow));
figure; plot(rois(:,temp) - 0.2*repmat(1:numToShow,t,1),'Color','k');
%% find peaks
if ~baselineexists
   rois = rois(601:end,:);
else
   rois = rois(301:600,:);
end
[t,n] = size(rois);
rois_mean = mean(rois,1);
if ~baselineexists
    rois_std = std(rois,[],1);
else
    rois_std = blStruct.stds;
end
rois_std(rois_std >  0.1) =  0.1;
rois_std(rois_std <  prctile(rois_std,5)) =  prctile(rois_std,5);
events = struct();
timeMinutes = size(rois,1)/60;

binarywhole = [];
for i=1:n
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


%%stats for things
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
[fp,name,~] = fileparts([dname fn]);

if ~baselineexists
    save([fp '\' name '_SGNstruct2.mat'],'SGNstruct');
else
    save([fp '\' name '_SGNstruct2' drug '.mat'],'SGNstruct');
end

    
%%
if ~baselineexists
    imgfilt = Kalman_Stack_Filter(single(img(:,:,601:end)));
else
    imgfilt = Kalman_Stack_Filter(single(img(:,:,301:600)));
end
    
generateActivityMovie(imgfilt, binarywhole,positiveIndices,[fp '\' name '_ActiveMovie'],[0 2500])
