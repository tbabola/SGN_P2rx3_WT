function [groupStruct] = groupedActivity(SGNstruct, threshold)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
    groupStruct = struct();
    for i = 1:size(SGNstruct,2)
        temp = reshape([SGNstruct(i).events.binary],size(SGNstruct(i).events(1).binary,2),[]);
        temp = sum(temp',1);
             
        [pks,locs,w] = findpeaks(temp,1:size(temp,2),'MinPeakProminence', threshold);
        %figure;
      %  findpeaks(temp,1:size(temp,2),'MinPeakProminence', threshold)
        groupStruct(i).freq = size(pks,2)/(size(temp,2)/60);
        groupStruct(i).meanROIs = mean(pks);
        
        count = 1; freqprc = [];
        for j = 80:5:100
            freqprc(count) = prctile([SGNstruct(i).events(~[SGNstruct(i).events.isEmpty]).frequency],j);
            count = count + 1;
        end
        groupStruct(i).freqPrcActive = freqprc;
    end
end

