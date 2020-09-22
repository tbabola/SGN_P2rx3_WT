function pairedFileList = findPairs(dir, key)
%UNTITLED2 This function returns pairs of matched files
%   Detailed explanation goes here
    if nargin < 2
        key = '';
    end
    fileList = loadFileList(dir);
    pairedFileList = {};
    
    finished = 0; count = 1; 
    
    while ~finished

        if ~isempty(fileList{count}) 
            [fdir] = fileparts(fileList{count});
            temp = regexp(fileList{count},'cochlea_\d+_base_','match');
            temp = [fdir '\' temp{1}];
            
            matchedFiles = fileList(contains(fileList,temp))
            if size(matchedFiles,1) > 1 && sum(contains(matchedFiles,key)) > 0
                pairedFileList = [pairedFileList; fileList(contains(fileList,temp))]
            end
            
            indices = find(contains(fileList,temp));
            for j = 1:size(indices,1)
                fileList(indices(j)) = {''};
            end
        end
        count = count + 1;
        if count > size(fileList,1)
            finished = 1;
        end
    end
end

