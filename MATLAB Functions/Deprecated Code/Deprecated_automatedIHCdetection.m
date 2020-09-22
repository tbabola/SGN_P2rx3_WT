clear all; close all;
addpath(genpath('..\MATLAB Functions'));
[fn, dname] = uigetfile('*.czi'); %prompt for file to analyze
%open image

%%
imgdata = bfopen([dname fn]);
imgdata = imgdata{1,1}; %this is where the imagedata is stored
[fp,name,~] = fileparts([dname fn]); %store filename and path for saving later
oimg = []; index = 1;
for i = 1:1:size(imgdata,1) %because there are two channels in your image, only select GCaMP signal
    oimg(:,:,index) = imgdata{i,1};
    index = index + 1;
end
oimg = bleachCorrect(oimg,2);
%%
meanImg = mean(oimg,3);
meanImg = (meanImg - mean(meanImg,'all'))/std(meanImg,[],'all');

figure; imagesc(meanImg);
fprintf('\n Select IHC search area using the mouse...\n  ')
roi = drawpolygon;
IHCmask = createMask(roi);

%%
%thresholding method
figure(1); 
maskedImg = meanImg.*IHCmask;
maskedImg(maskedImg == 0) = 1000;
imagesc(imbinarize(imcomplement(maskedImg),'adaptive','Sensitivity',1));
%imagesc(imbinarize(imcomplement(meanImg.*IHCmask),'adaptive','Sensitivity',1));
binarized = imbinarize(imcomplement(meanImg.*IHCmask),'adaptive','Sensitivity',1);
structs = regionprops(binarized,'Centroid','Area','PixelList');
structs = structs([structs.Area] < 200 & [structs.Area] > 5)
roiHandles = {};
figure(2); imagesc(meanImg);
for i=1:size(structs,1)
    roiHandles{i} = drawpoint('Position',structs(i).Centroid);
end
% roi = images.roi.Rectangle(gca,'Position',[12 48 20 20],'FixedAspectRatio',true);
% templateMask = find(createMask(roi));

disp('Add ROIs manually, then hit escape.')
points = []; roi = []; addedPoints = 0;
while 1
    temp = drawpoint;
    if isempty(temp.Position)
        break
    else
        roiHandles{end+1} = temp;
        addedPoints = 1;
    end
end

for i=1:size(roiHandles,2)
    points(i,:) = roiHandles{i}.Position;
end

if addedPoints
    tempStr = input('Sort by x or y?','s');
    if strcmp(tempStr,'x')
        [~,temp] = sort(points(:,1));
        points = points(temp,:);
    elseif strcmp(tempStr,'y')
        [~,temp] = sort(points(:,2));
        points = points(temp,:);
    end
end
 
oPoints = points;

%%
f = fit(points(:,1),points(:,2),'poly4');
tempStr = [min(points(:,1)):0.5:max(points(:,1))]';
polypoints = [tempStr f(tempStr)];

roi = drawpolyline('Position',polypoints(1:15:end,:));
input('Drag line to the medial side of the IHCs, then press enter.')
newpolypoints = polypoints - (polypoints(16,:) - roi.Position(2,:));
%%
figure(2); imagesc(meanImg);
bottomPos = [];
distTowardsLine = 7; %number of pixels to move towards the line

for i = 1:size(points,1)
       %find closest Pixel
%     tempList = structs(i).PixelList;
%     dist = pdist2(tempList,polypoints);
%     [~,I]=min(min(dist,[],2));
%     bottomPos(i,:) = tempList(I,:);
%     drawpoint('Position',bottomPos(i,:),'Color','g');

       %go from centroid towards drawn line
      tempCentroid = points(i,:);
      dist = pdist2(tempCentroid,newpolypoints);
      [~,I]= min(dist);
      tempPoint = newpolypoints(I,:);
      tempslope = tempCentroid - tempPoint;
      %a little geometry
      tempSq = tempslope.^2;
      tempSqAdd = sum(tempSq);
      unitToMove = sqrt(distTowardsLine^2/tempSqAdd);
      tempPoint = tempCentroid - tempslope.*unitToMove;
      bottomPos(i,:) = tempPoint;
%     drawpoint('Position',bottomPos(i,:),'Color','g');
end

figure(2); imagesc(meanImg);
for i = 1:size(bottomPos,1)
    drawpoint('Position',points(i,:),'Color','r');
    drawpoint('Position',bottomPos(i,:),'Color','g');
end

%%
% movAvg = 5;
% slope = 1;
% points = bottomPos;
% points(:,2) = 512 - points(:,2) ;
% 
% for i=1:size(points,1)
%     startPt = i - floor(movAvg/2);
%     endPt = i + floor(movAvg/2);
%     startPt;
%     if startPt < 1
%         startPt = 1;
%     elseif endPt > size(points,1)
%         endPt = size(points,1);
%     end
%     p = polyfit(points(startPt:endPt,1),points(startPt:endPt,2),1);
%     slope = [slope; p(1)];
% end
% 
% angles = atan(slope) * 180/pi;
% angles = angles + 90;

%%
roisignal = zeros(size(oimg,3),size(points,1));
I = {}; roihandles = {};
centers = points;
points = bottomPos;
figure(2); imagesc(meanImg);
for i=1:size(points,1)
    roihandles{i} = drawellipse('Center',points(i,:),'SemiAxes',[3 3]);
end

input('Adjust any points and then hit enter at the command line.');

for i=1:size(points,1)
    h = roihandles{i};
    finalPoints(i,:) = h.Position;
    mask = h.createMask;
    I{i} = find(mask);
end

for i=1:size(oimg,3)
    temp = oimg(:,:,i);
    for j=1:size(I,2)
        tempI = I{j};
        roisignal(i,j) = mean(temp(tempI),'all');
    end 
end
truesize; 

%%
figure(3); 
% for i = 1:10
%     randomIdx = round(rand(1)*size(roisignal,2));
%     if randomIdx == 0
%         randomIdx = 1;
%     end
%     hold off; plot(smooth(roisignal(:,i))); hold on;
%     %meanVal = mean(roisignal(:,randomIdx));
%     prcVal = prctile(smooth(roisignal(:,randomIdx)),5);
%     modeVal = mode(roisignal(:,randomIdx));
%     %plot([0 600],[meanVal meanVal]);
%     plot([0 600],[prcVal prcVal]);
%     pause(1);
% end

prcF = prctile(roisignal,5);
roisignalNorm = (roisignal - prcF) ./prcF;
plotOffset(roisignalNorm, 0.2);
for i = 1:size(roisignalNorm,2)
    roisignalNorm(:,i) = smooth(roisignalNorm(:,i),5);
end

%%
medians = []; stds = []; window = 15;
for i = 1:size(roisignalNorm,1)-window
   medians(i,:) =  median(roisignalNorm(i:i+window-1,:));
   stds(i,:) = std(roisignalNorm(i:i+window-1,:),[],1);
end

finMedian = median(medians);
finStds = median(stds);
thr = finMedian + 4*finStds;

binary = roisignalNorm > thr;
generateActivityMovieIHCs(Kalman_Stack_Filter(single(oimg)), binary, centers,'testOfIHCs',[min(oimg,[],'all') 20000]);

%%
%%enhance
meanImg = mean(oimg,3);
meanImg = (meanImg - mean(meanImg,'all'))/std(meanImg,[],'all');
figure; imagesc(meanImg);

while 1
    roi = drawpoint();
    if isempty(roi.Position)
        break
    end
    pos = roi.Position;
    val = meanImg(round(pos(1)),round(pos(2)))
    tempImg = meanImg;
    tempImg(tempImg < val) = tempImg(tempImg < val) * 0.9;
    imagesc(tempImg);
    meanImg = tempImg;
end



% %%
% meanImgMask = meanImg .*IHCmask;
% template = meanImg(templateMask);
% template = reshape(template,20,20);
% bigTemp = NaN(30,30);
% bigTemp(6:25,6:25) = template;
% %%
% test=[];
% count = 1;
% for i= 0:5:90
%     angle = i;
%     tempOut = imrotate(bigTemp,-angle,'bilinear','crop');
%     %imagesc(tempOut);
%     
%     test(:,:,count) = normxcorr2(tempOut,meanImgMask);
%     imagesc(normxcorr2(tempOut,meanImgMask));
%     drawnow;
%     count =count + 1;
% end
% 
% 
% %%okay
% % sizeOfPattern = 13;
% % pattern = zeros(sizeOfPattern,sizeOfPattern);
% % patternCenter = ceil(sizeOfPattern/2);
% % pattern(patternCenter,patternCenter) = 1;
% % out = bwdist(pattern) <= 3.5;
% % out = imresize(out,[20 10]);
% % figure(2); imagesc(out);
% % 
% %  figure(3); 
% %  count = 1;
% % for i= 0:10:90
% %     angle = i;
% %     tempOut = imrotate(out,-angle,'bilinear','crop');
% %     tempOut = tempOut == 0;
% %     %imagesc(tempOut);
% %     
% %     test(:,:,count) = conv2(meanImg,tempOut,'same');
% %     imagesc(xcorr2(meanImg,tempOut));
% %     drawnow;
% %     count =count + 1;
% end
%figure(3); imagesc(conv2(meanImg,out))
