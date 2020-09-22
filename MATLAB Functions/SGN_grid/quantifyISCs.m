%%
%playing around
imgx = 512;
imgy = 512;
img = zeros(imgx,imgy);
img = imrotate(img,45);

[m,n] = size(img);
grid = getGrid(m,n,10);
mask = ISCstruct.ISCmask;
posGrid = getPositiveGrid(grid,mask.Position);

%%
gridIndex = [];
count = 1;
for i = 1:ceil(m/10)
    for j = 1:ceil(n/10)
        gridIndex(count,:) = [i j];
        count = count + 1;
    end
end

[~,idx] = ismember(posGrid,grid,'rows')

%%
scatPoints = [];

for i = 1:size(idx)
    hits = find(ISCstruct.events(i).binary);
    scatPoints = [scatPoints; repmat(gridIndex(idx(i),1),size(hits,2),1) repmat(gridIndex(idx(i),2),size(hits,2),1) hits'];
end
%%
figure; scatter3(scatPoints(:,1),scatPoints(:,2),scatPoints(:,3));

%%
T = clusterdata(scatPoints(1:1250,:),'MaxClust',40);
X = scatPoints(1:1250,:);
scatter3(X(:,1),X(:,2),X(:,3),100,T,'filled')
title('Result of Clustering');

%%
Z = linkage(X);
dendrogram(Z)

%%
%try a slightly different approach
roiSignal = ISCstruct.rois;
medSig = median(roiSignal);
stdSig = std(roiSignal,1);

boohiss = roiSignal > medSig + 3*stdSig;

%%
imgbinary = zeros(ceil(size(img,1)/10),ceil(size(img,2)/10),1250);

for i = 1:size(idx)
    imgbinary(gridIndex(idx(i),1),gridIndex(idx(i),2),:) = boohiss(:,i);
end

imgbinary = imgbinary > 0;
labels = bwlabeln(imgbinary,26);
count = 0; number = []; wss = [];

newlabels = zeros(size(labels));
labelcount = 1;
for i=1:max(labels,[],'all') 
    [i1,i2,i3] = ind2sub(size(labels),find(labels(:)==i)); %get x,y,t coordinates for label i
    
    if max(i3)-min(i3) <= 1 | size(i1,1) <= 4 %scrub any events that last less than 3 frames or has less than 4 connected squares
        count = count + 1;
    elseif size(i1,1) > 500
        temp = (labels == i) .* labels;
        singleLabel = labelsToRois(temp, gridIndex,idx,size(boohiss));
        generateActivityMovie_labels2(img(:,:,min(i3):max(i3)),singleLabel(min(i3):max(i3),:)',posGrid,[],[475 16000]);
        drawnow;
        
       x = queryEventNum();
       
       satisfied = 0;
       
       xscale = 1;
       yscale = 1;
       tscale = 2;
       while ~satisfied
           if x == 1
               satisfied = 1;
               newlabels(labels == i) = labelcount;
               labelcount = labelcount + 1;
           else
               indices = kmeans([i1/xscale i2/yscale i3/tscale],x); %cluster data into x number of groups
               temp(sub2ind(size(labels),i1,i2,i3)) = indices; %store these clusters temporarily
               singleLabel = labelsToRois(temp, gridIndex,idx,size(boohiss)); %convert into ROIs for generate activity movie function
               generateActivityMovie_labels2(img(:,:,min(i3):max(i3)),singleLabel(min(i3):max(i3),:)',posGrid,[],[475 16000]);
               
               satisfied = querySatisfaction();
               
               if ~satisfied
                   [xscale, yscale, tscale] = queryScale(xscale,yscale,tscale);
               else
                   for j = 1:max(indices)
                        newlabels(temp == j) = labelcount;
                        labelcount = labelcount + 1;
                   end
               end
           end
       end
    else %if small, just add label into new labels
        newlabels(labels == i) = labelcount;
        labelcount = labelcount + 1;
    end
end
disp([num2str(count) ' events were scrubbed.'])
%%
labelRoi = labelsToRois(newlabels, gridIndex, idx, size(boohiss))
sum
%%
generateActivityMovie_labels(Kalman_Stack_Filter(single(img)),labelRoi',posGrid,'test',[475 16000]);

%% analysis
ISCstruct.activeArea = sum(sum(labelRoi,1) > 1) / size(labelRoi,2)
for i = 1:max(newlabels,[],'all')
    [i1,i2,i3] = ind2sub(size(newlabels),find(newlabels(:)==i)); %get x,y,t coordinates for label i
    ISCstruct.event(i).timeStart = min(i3);
    ISCstruct.event(i).timeEnd = max(i3);
    ISCstruct.event(i).eventDuration = max(i3) - min(i3);
    ISCstruct.event(i).area = size(unique([i1 i2],'rows'),1)
end


%%
function labelRoi = labelsToRois(labels, gridIndex, posIndex, sizeRois)
    labelRoi = zeros(sizeRois);
    for i = 1:size(posIndex)
        labelRoi(:,i) = squeeze(labels(gridIndex(posIndex(i),1),gridIndex(posIndex(i),2),:));
    end
end

function numEvents = queryEventNum()
    x = input('How many events do you see? ');
    
       while 1
           if isnumeric(x) && x > 0
               numEvents = x;
               return
           else
             x = input('How many events do you see? (Enter a postive number) ');
           end
       end
end

function satisfied = querySatisfaction()
    satisfied = input('Satisfied?');
    
    while 1
       if isnumeric(satisfied) && satisfied >= 0 && satisfied <=1
           return
       else
         x = input('Are you satisfied? (0) No (1) Yes');
       end
    end
end

function [xscale,yscale,tscale] = queryScale(xscale,yscale,tscale)
    xscale = xscale; yscale = yscale; tscale = tscale;
    x = input(['Which scale would you like to alter to adjust clustering?\n' ... 
        '(1) x, Current value: ' num2str(xscale) '\n' ...
        '(2) y, Current value: ' num2str(yscale) '\n' ...
        '(3) t, Current value: ' num2str(tscale) '\n']);
    
    while 1
       if x > 0 && x <= 3
           if x == 1
                xscale = input('What is the new scale? ');
                while 1
                     if isnumeric(xscale)
                           return
                     else
                        xscale = input('What is the new scale? (Must be numeric) ');
                     end
                end
           elseif x == 2
                yscale = input('What is the new scale?');
                while 1
                     if isnumeric(yscale)
                           return
                     else
                        yscale = input('What is the new scale? (Must be numeric) ');
                     end
                end
           elseif x == 3
                tscale = input('What is the new scale?');
                while 1
                     if isnumeric(tscale)
                           return
                     else
                        tscale = input('What is the new scale? (Must be numeric) ');
                     end
                end
           end
       else
         x = input(['Which scale would you like to alter to adjust clustering?\n' ... 
        '(1) x, Current value: ' num2str(xscale) '\n' ...
        '(2) y, Current value: ' num2str(yscale) '\n' ...
        '(3) t, Current value: ' num2str(tscale) '\n']);
       end
   end
end