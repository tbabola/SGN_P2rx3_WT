function tempCells = loadCellStructs(dname)
    
    if class(dname) == 'char'
        fileList = loadFileList(dname);
    elseif class(dname) == 'cell'
        fileList = dname;
    end
    
    for i=1:size(fileList,1)
        s = load(fileList{i});
        fns = fieldnames(s);
        tempCells(i)= s.(fns{1});
    end
end
