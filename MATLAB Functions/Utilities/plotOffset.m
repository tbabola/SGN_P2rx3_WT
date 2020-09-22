function [h] = plotOffset(data,offset)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    [t,n] = size(data); %t in first dimension
    h = figure; plot(data - offset*repmat(1:n,t,1),'Color','k');
end

